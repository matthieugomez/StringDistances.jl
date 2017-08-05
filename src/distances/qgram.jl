##############################################################################
##
## Define a type that iterates through q-grams of a string
##
##############################################################################

struct QGramIterator{S <: AbstractString, T <: Integer}
	s::S # grapheme
	l::Int # length of string
	q::T # length of q-grams
end

function Base.start(qgram::QGramIterator)
	(1, qgram.l < qgram.q ? endof(qgram.s) + 1 : chr2ind(qgram.s, qgram.q))
end
function Base.next(qgram::QGramIterator, state)
	istart, iend = state
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
function Base.done(qgram::QGramIterator, state)
	istart, idend = state
	done(qgram.s, idend)
end
Base.eltype(qgram::QGramIterator) = SubString{typeof(qgram.s)}
Base.length(qgram::QGramIterator) = max(qgram.l - qgram.q + 1, 0)
function Base.collect(qgram::QGramIterator)
	x = Array{eltype(qgram)}(length(qgram))
	i = 0
	for q in qgram
		i += 1
		@inbounds x[i] = q
	end
	x
end
Base.sort(qgram::QGramIterator) = sort!(collect(qgram))

##############################################################################
##
## For each element in union{v1, v2}, this iterator output numbers of times it appears in v1 and the number of times it appears in v2
## v1 and v2 must be sorted vectors
##
##############################################################################

struct CountInterator{T1 <: AbstractVector, T2 <: AbstractVector}
	v1::T1
	v2::T2
end
Base.start(s::CountInterator) = (1, 1)

function Base.next(s::CountInterator, state)
	state1, state2 = state
	iter1 = done(s.v2, state2)
	iter2 = done(s.v1, state1)
	if iter1
		@inbounds x1 = s.v1[state1]
	elseif iter2
		@inbounds x2 = s.v2[state2]
	else
		@inbounds x1 = s.v1[state1]
		@inbounds x2 = s.v2[state2]
		iter1 = x1 <= x2
		iter2 = x2 <= x1
	end
	nextstate1 = iter1 ? searchsortedlast(s.v1, x1, state1, length(s.v1), Base.Forward) + 1 : state1
	nextstate2 = iter2 ? searchsortedlast(s.v2, x2, state2, length(s.v2), Base.Forward) + 1 : state2
	((nextstate1 - state1, nextstate2 - state2), (nextstate1, nextstate2))
end

function Base.done(s::CountInterator, state) 
	state1, state2 = state
	done(s.v2, state2) && done(s.v1, state1)
end

##############################################################################
##
## Distance on strings is computed by set distance on qgram sets
##
##############################################################################
abstract type AbstractQGram <: SemiMetric end

function evaluate(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
	sort1 = sort(QGramIterator(s1, length(s1), dist.q))
	sort2 = sort(QGramIterator(s2, length(s2), dist.q))
	evaluate(dist, CountInterator(sort1, sort2))
end

##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

struct QGram{T <: Integer} <: AbstractQGram
	q::T
end
QGram() = QGram(2)

function evaluate(dist::QGram, countiterator)
	n = 0
	for (n1, n2) in countiterator
		n += abs(n1 - n2)
	end
	n
end

##############################################################################
##
## cosine 
##
## 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
##############################################################################

struct Cosine{T <: Integer} <: AbstractQGram
	q::T
end
Cosine() = Cosine(2)

function evaluate(dist::Cosine, countiterator)
	norm1, norm2, prodnorm = 0, 0, 0
	for (n1, n2) in countiterator
		norm1 += n1^2
		norm2 += n2^2
		prodnorm += n1 * n2
	end
	1.0 - prodnorm / (sqrt(norm1) * sqrt(norm2))
end

##############################################################################
##
## Jaccard
##
## Denote Q(s, q) the set of tuple of length q in s
## 1 - |intersect(Q(s1, q), Q(s2, q))| / |union(Q(s1, q), Q(s2, q))|
##
##############################################################################

struct Jaccard{T <: Integer} <: AbstractQGram
	q::T
end
Jaccard() = Jaccard(2)

function evaluate(dist::Jaccard, countiterator)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in countiterator
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / (ndistinct1 + ndistinct2 - nintersect)
end

##############################################################################
##
## SorensenDice
##
## 1 - 2 * |intersect(Q(s1, q), Q(s2, q))| / (|Q(s1, q)| + |Q(s2, q))|)
##############################################################################

struct SorensenDice{T <: Integer} <: AbstractQGram
	q::T
end
SorensenDice() = SorensenDice(2)

function evaluate(dist::SorensenDice, countiterator)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in countiterator
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - 2.0 * nintersect / (ndistinct1 + ndistinct2)
end

##############################################################################
##
## overlap
##
## 1 -  |intersect(Q(s1, q), Q(s2, q))| / min(|Q(s1, q)|, |Q(s2, q)))
##############################################################################

struct Overlap{T <: Integer} <: AbstractQGram
	q::T
end
Overlap() = Overlap(2)

function evaluate(dist::Overlap, countiterator)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in countiterator
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / min(ndistinct1, ndistinct2)
end

