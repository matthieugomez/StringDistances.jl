
##############################################################################
##
## Define a type that iterates through q-grams of a string
##
##############################################################################

struct QGramIterator{S <: AbstractString, N}
	s::S # grapheme
	l::Int # length of string
end
param(x::QGramIterator{S, N}) where {S, N} = N

function Base.iterate(qgram::QGramIterator{S, N}, 
	state = (1, qgram.l < N ? lastindex(qgram.s) + 1 : nextind(qgram.s, 0, N))) where {S, N}
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = qgram.s[istart:iend]
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.length(qgram::QGramIterator{S, N}) where {S, N} = max(qgram.l - N + 1, 0)

Base.eltype(qgram::QGramIterator) = String

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

# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
function CountIteratorDictionary(s1::QGramIterator{S1, N}, s2::QGramIterator{S2, N}) where {S1, S2, N}
	K = String
	d = Dict{K, NTuple{2, Int}}()
	sizehint!(d, length(s1) + length(s2))
	for ch10 in s1
		ch1 = convert(K, ch10)
		!isequal(ch1, ch10) && throw(ArgumentError("$(limitrepr(ch10)) is not a valid key for type $K"))
		index = Base.ht_keyindex2!(d, ch1)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = ch1
			@inbounds d.vals[index] = (d.vals[index][1] + 1, 0)
		else
			Base._setindex!(d, (1, 0), ch1, -index)
		end
	end
	for ch20 in s2
		ch2 = convert(K, ch20)
		!isequal(ch2, ch20) && throw(ArgumentError("$(limitrepr(ch20)) is not a valid key for type $K"))
		index = Base.ht_keyindex2!(d, ch2)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = ch2
			@inbounds d.vals[index] = (d.vals[index][1], d.vals[index][2] + 1)
		else
			Base._setindex!(d, (0, 1), ch2, -index)
		end
	end
	return values(d)
end



##############################################################################
##
## Distance on strings is computed by set distance on qgram sets
##
##############################################################################
abstract type AbstractQGram{N} <: SemiMetric end
param(x::AbstractQGram{N}) where N = N


function evaluate(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
	evaluate(dist, 
		CountIteratorDictionary(QGramIterator{typeof(s1), param(dist)}(s1, length(s1)), 
		QGramIterator{typeof(s2), param(dist)}(s2, length(s2))))
end

##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

struct QGram{N} <: AbstractQGram{N} end

QGram(x::Integer) = QGram{x}()

function evaluate(dist::QGram, countiterator)
	n = 0
	for (n1, n2) in countiterator
		n += abs(Int(n1) - Int(n2))
	end
	n
end

##############################################################################
##
## cosine 
##
## 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
##############################################################################

struct Cosine{N} <: AbstractQGram{N} end

Cosine(n::Integer = 2) = Cosine{n}()

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

struct Jaccard{N} <: AbstractQGram{N} end

Jaccard(n::Integer = 2) = Jaccard{n}()

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

struct SorensenDice{N} <: AbstractQGram{N} end

SorensenDice(n::Integer = 2) = SorensenDice{n}()

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

struct Overlap{N} <: AbstractQGram{N} end

Overlap(n::Integer = 2) = Overlap{n}()

function evaluate(dist::Overlap, countiterator)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in countiterator
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / min(ndistinct1, ndistinct2)
end

