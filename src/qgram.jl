
##############################################################################
##
## Define a type that iterates through q-grams of a string
##
############################################################################
struct QGramIterator{S <: AbstractString}
	s::S # grapheme
	l::Int # length of string
	q::Int # Length of Qgram
end

function Base.iterate(qgram::QGramIterator, 
	state = (1, qgram.l < qgram.q ? ncodeunits(qgram.s) + 1 : nextind(qgram.s, 0, qgram.q)))
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.length(qgram::QGramIterator) = max(qgram.l - qgram.q + 1, 0)
Base.eltype(qgram::QGramIterator{SubString{S}}) where {S} = SubString{S}
Base.eltype(qgram::QGramIterator{S}) where {S} = SubString{S}

"""
Return an iterator that iterates on the QGram of the string

### Arguments
* `s::AbstractString`
* `q::Integer`: length of qgram

## Examples
```julia
using StringDistances
for x in qgram("hello", 2)
	@show x
end
```
"""
function qgram(s::AbstractString, q::Integer)
	QGramIterator{typeof(s)}(s, length(s), q)
end

##############################################################################
##
## 
## 
## 
##
##############################################################################
"""
   count_map(x1, x2)

For two iterators `x1` and `x2`, `count_map(x1, x2)`  returns an dictionary 
that returns,  for each element in `x1` or `x2`, a tuple with the numbers of 
times it appears in `x1` and the number of times it appears in `x2`
"""
function count_map(s1, s2)
	K = promote_type(eltype(s1), eltype(s2))
	d = Dict{K, Tuple{Int, Int}}()
	# I use a faster way to change a dictionary key
	# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
	sizehint!(d, length(s1) + length(s2))
	for x1 in s1
		index = Base.ht_keyindex2!(d, x1)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x1
			@inbounds d.vals[index] = (d.vals[index][1] + 1, 0)
		else
			@inbounds Base._setindex!(d, (1, 0), x1, -index)
		end
	end
	for x2 in s2
		index = Base.ht_keyindex2!(d, x2)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x2
			@inbounds d.vals[index] = (d.vals[index][1], d.vals[index][2] + 1)
		else
			@inbounds Base._setindex!(d, (0, 1), x2, -index)
		end
	end
	return d
end

#= Trie
function count_map(s1, s2)
	d = Trie{Tuple{Int, Int}}()
	for ch1 in s1
		node = d
		for char in ch1
			if !haskey(node.children, char)
	            node.children[char] = Trie{Tuple{Int, Int}}()
	        end
	        node = node.children[char]
	    end
	    node.value = node.is_key ? (node.value[1] + 1, 0) : (1, 0)
	    node.is_key = true
	end
	for ch2 in s2
		node = d
		for char in ch2
			if !haskey(node.children, char)
	            node.children[char] = Trie{Tuple{Int, Int}}()
	        end
	        node = node.children[char]
	    end
	    node.value = node.is_key ? (node.value[1], node.value[2]+ 1) : (0, 1)
	    node.is_key = true
	end
	return iterator(d)
end
function iterator(t::Trie, found = Tuple{Int, Int}[])
    if t.is_key
    	t.is_key = false
    	push!(found, t.value)
    else
    	for k in values(t.children)
    		iterator(k, found) 
    	end
    end
    return found
end
=#


##############################################################################
##
## Distance on strings is computed by set distance on qgram sets
##
##############################################################################
abstract type AbstractQGramDistance <: SemiMetric end

function evaluate(dist::AbstractQGramDistance, s1::AbstractString, s2::AbstractString)
	x = count_map(qgram(s1, dist.q), qgram(s2, dist.q))
	evaluate(dist, x)
end

##############################################################################
##
## q-gram 
##
##############################################################################
"""
For an AbstractString s, denote v(s) the vector on the space of q-grams of length N, that contains the number of times a q-gram appears in s
The q-gram distance is ||v(s1) - v(s2)||
"""

"""
	QGram(q::Int)

Creates a QGram metric.

The distance corresponds to

``||v(s1, q) - v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, that contains the number of times a q-gram appears for the string s
"""
struct QGram <: AbstractQGramDistance
	q::Int
end

function evaluate(dist::QGram, count_dict)
	n = 0
	for (n1, n2) in values(count_dict)
		n += abs(n1 - n2)
	end
	n
end

##############################################################################
##
## cosine 
##
## 
##############################################################################
"""
	Cosine(q::Int)

Creates a Cosine metric.

The distance corresponds to

`` 1 - v(s1, q).v(s2, q)  / ||v(s1, q)|| * ||v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, that contains the number of times a q-gram appears for the string s
"""
struct Cosine <: AbstractQGramDistance
	q::Int
end

function evaluate(dist::Cosine, count_dict)
	norm1, norm2, prodnorm = 0, 0, 0
	for (n1, n2) in values(count_dict)
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
##############################################################################
"""
	Jaccard(q::Int)

Creates a Jaccard metric.

The distance corresponds to 

``1 - |Q(s1, q) ∩ Q(s2, q)| / |Q(s1, q) ∪ Q(s2, q))|``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Jaccard <: AbstractQGramDistance
	q::Int
end

function evaluate(dist::Jaccard, count_dict)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in values(count_dict)
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
##############################################################################
"""
	SorensenDice(q::Int)

Creates a SorensenDice metric

The distance corresponds to  

``1 - 2 * |Q(s1, q) ∩ Q(s2, q)|  / (|Q(s1, q)| + |Q(s2, q))|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct SorensenDice <: AbstractQGramDistance
	q::Int
end

function evaluate(dist::SorensenDice, count_dict)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in values(count_dict)
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
"""
	Overlap(q::Int)

Creates a Overlap metric

The distance corresponds to  

``1 - |Q(s1, q) ∩ Q(s2, q)|  / min(|Q(s1, q)|, |Q(s2, q)|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Overlap <: AbstractQGramDistance
	q::Int
end

function evaluate(dist::Overlap, count_dict)
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in values(count_dict)
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / min(ndistinct1, ndistinct2)
end

