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


function Base.iterate(qgram::QGramIterator, 
	state = (1, qgram.l < qgram.q ? lastindex(qgram.s) + 1 : nextind(qgram.s, 0, qgram.q)))
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.eltype(qgram::QGramIterator{S}) where {S} = SubString{S}
Base.eltype(qgram::QGramIterator{S}) where {S <: SubString} = S
Base.length(qgram::QGramIterator) = max(qgram.l - qgram.q + 1, 0)

##############################################################################
##
## CountedIterator that use Binary Search
##
## For each element in union{v1, v2}, this iterator output numbers of times it appears in v1 and the number of times it appears in v2
## v1 and v2 must be sorted vectors
##
##############################################################################
struct CountIteratorBinary{T1, T2}
	v1::Vector{T1}
	v2::Vector{T2}
end

function Base.collect(qgram::QGramIterator)
	x = Array{eltype(qgram)}(undef, length(qgram))
	i = 0
	for q in qgram
		i += 1
		@inbounds x[i] = q
	end
	x
end
Base.sort(qgram::QGramIterator) = sort!(collect(qgram))


function CountIteratorBinary(s1::QGramIterator, s2::QGramIterator)
	CountIteratorBinary(sort(s1), sort(s2))
end

function Base.iterate(s::CountIteratorBinary, state = (1, 1))
	state1, state2 = state
	iter1 = state2 > length(s.v2)
	iter2 = state1 > length(s.v1)
	iter2 && iter1 && return nothing
	if iter1
		x1 = s.v1[state1]
	elseif iter2
		x2 = s.v2[state2]
	else
		x1 = s.v1[state1]
		x2 = s.v2[state2]
		iter1 = x1 <= x2
		iter2 = x2 <= x1
	end
	nextstate1 = iter1 ? searchsortedlast(s.v1, x1, state1, length(s.v1), Base.Forward) + 1 : state1
	nextstate2 = iter2 ? searchsortedlast(s.v2, x2, state2, length(s.v2), Base.Forward) + 1 : state2
	((nextstate1 - state1, nextstate2 - state2), (nextstate1, nextstate2))
end


##############################################################################
##
## CountedIterator that use Dictionary
##
## For each element in union{v1, v2}, this iterator output numbers of times it appears in v1 and the number of times it appears in v2
## v1 and v2 must be sorted vectors
##
##############################################################################
struct CountIteratorDictionary{T}
	d::T
end

function CountIteratorDictionary(s1::QGramIterator, s2::QGramIterator)
	d = Dict{eltype(s1), Tuple{Int, Int}}()
	for ch1 in s1
		if haskey(d, ch1)
			t = d[ch1]
			d[ch1] = (t[1] + 1, 0)
		else
			d[ch1] = (1, 0)
		end
	end
	for ch2 in s2
		if haskey(d, ch2)
			t = d[ch2]
			d[ch2] = (t[1], t[2] + 1)
		else
			d[ch2] = (0, 1)
		end
	end
	return values(d)
end



##############################################################################
##
## Distance on strings is computed by set distance on qgram sets
##
##############################################################################
abstract type AbstractQGram <: SemiMetric end

function evaluate(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
	evaluate(dist, 
		CountIteratorBinary(QGramIterator(s1, length(s1), dist.q), 
		QGramIterator(s2, length(s2), dist.q)))
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

