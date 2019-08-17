
##############################################################################
##
## Define a type that iterates through q-grams of a string
##
##############################################################################
# N is the number of characters in the QGram
struct QGramIterator{S <: AbstractString}
	s::S # grapheme
	l::Int # length of string
	N::Int # Length of Qgram
end

"""
Return an iterator that iterates on the QGram of the string

### Arguments
* `s::AbstractString`
* `n::Int`: length of qgram

## Examples
```julia
using StringDistances
for x in qgram_iterator("hello", 2)
	@show x
end
```
"""
qgram_iterator(s::AbstractString, n::Int) = QGramIterator{typeof(s)}(s, length(s), n)

function Base.iterate(qgram::QGramIterator, 
	state = (1, qgram.l < qgram.N ? ncodeunits(qgram.s) + 1 : nextind(qgram.s, 0, qgram.N)))
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = qgram.s[istart:iend]
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.length(qgram::QGramIterator) = max(qgram.l - qgram.N + 1, 0)
Base.eltype(qgram::QGramIterator{SubString{S}}) where {S} = S
Base.eltype(qgram::QGramIterator{S}) where {S} = S

##############################################################################
##
## For two iterators x1 x2, count_map(x1, x2) returns an iterator 
## that returns,  for each element in union{x1, x2}, the numbers of 
## times it appears in x1 and the number of times it appears in x2
##
##############################################################################

# I use a faster way to change a dictionary key
# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
function count_map(s1, s2)
	K = Union{eltype(s1), eltype(s2)}
	d = Dict{K, NTuple{2, Int}}()
	sizehint!(d, length(s1) + length(s2))
	for ch1 in s1
		index = Base.ht_keyindex2!(d, ch1)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = ch1
			@inbounds d.vals[index] = (d.vals[index][1] + 1, 0)
		else
			@inbounds Base._setindex!(d, (1, 0), ch1, -index)
		end
	end
	for ch2 in s2
		index = Base.ht_keyindex2!(d, ch2)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = ch2
			@inbounds d.vals[index] = (d.vals[index][1], d.vals[index][2] + 1)
		else
			@inbounds Base._setindex!(d, (0, 1), ch2, -index)
		end
	end
	return values(d)
end

##############################################################################
##
## Distance on strings is computed by set distance on qgram sets
##
##############################################################################
abstract type AbstractQGramDistance <: SemiMetric end

function evaluate(dist::AbstractQGramDistance, s1::AbstractString, s2::AbstractString)
	evaluate(dist, 
		count_map(qgram_iterator(s1, dist.N), qgram_iterator(s2, dist.N)))
end

##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

struct QGram <: AbstractQGramDistance
	N::Int
end

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

struct Cosine <: AbstractQGramDistance
	N::Int
end

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

struct Jaccard <: AbstractQGramDistance
	N::Int
end

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

struct SorensenDice <: AbstractQGramDistance
	N::Int
end

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

struct Overlap <: AbstractQGramDistance
	N::Int
end

function evaluate(dist::Overlap, countiterator)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in countiterator
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / min(ndistinct1, ndistinct2)
end

