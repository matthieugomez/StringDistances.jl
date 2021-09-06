struct QGramIterator{S <: Union{AbstractString, AbstractVector}}
	s::S   # Collection
	q::Int # Length of Qgram
	function QGramIterator{S}(s, q) where {S <: Union{AbstractString, AbstractVector}}
		q > 0 || throw(ArgumentError("The qgram length must be higher than zero"))
		new(s, q)
	end
end
function QGramIterator(s::Union{AbstractString, AbstractVector}, q::Integer)
	QGramIterator{typeof(s)}(s, q)
end
Base.length(qgram::QGramIterator) = max(length(qgram.s) - qgram.q + 1, 0)

# q-grams of AbstractString
function Base.iterate(qgram::QGramIterator{<: AbstractString}, 
	state = (1, nextind(qgram.s, 0, qgram.q)))
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.eltype(qgram::QGramIterator{SubString{S}}) where {S} = SubString{S}
Base.eltype(qgram::QGramIterator{S}) where {S <: AbstractString} = SubString{S}
qgrams(s::AbstractString, q::Integer) = QGramIterator(s, q)


#q-grams of AbstractVector
# Alternatively, I could also use partition in IterTools but it creates a vector for each iteration
# so it does not seem to be worth it.
function Base.iterate(qgram::QGramIterator{<: AbstractVector}, state = firstindex(qgram.s))
	state + qgram.q - 1 > lastindex(qgram.s) && return nothing
	view(qgram.s, state:(state + qgram.q - 1)), state + 1
end
Base.eltype(qgram::QGramIterator{<: AbstractVector}) = typeof(first(qgram))
qgrams(s::AbstractVector, q::Integer) = QGramIterator(s, q)
qgrams(s, q::Integer) = QGramIterator(collect(s), q)

@doc """
Return an iterator corresponding to the the q-gram of an iterator. 
When the iterator is a String, qgrams are SubStrings.

### Arguments
* `s` iterator
* `q::Integer`: length of q-gram

## Examples
```julia
for x in qgrams("hello", 2)
	println(x)
end
```
""" 
qgrams

# For two iterators s1 and s2, that define a length and eltype method,
# this returns an iterator that,
# for each element in s1 ∪ s2, returns (numbers of times it appears in s1, numbers of times it appears in s2)
function _count(qgrams1, qgrams2)
	K = promote_type(eltype(qgrams1), eltype(qgrams2))
	d = Dict{K, Tuple{Int, Int}}()
	sizehint!(d, length(qgrams1) + length(qgrams2))
	# I use a faster way to change a dictionary key
	# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
	for x1 in qgrams1
		index = Base.ht_keyindex2!(d, x1)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x1
			@inbounds d.vals[index] = (d.vals[index][1] + 1, 0)
		else
			@inbounds Base._setindex!(d, (1, 0), x1, -index)
		end
	end
	for x2 in qgrams2
		index = Base.ht_keyindex2!(d, x2)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x2
			@inbounds d.vals[index] = (d.vals[index][1], d.vals[index][2] + 1)
		else
			@inbounds Base._setindex!(d, (0, 1), x2, -index)
		end
	end
	return values(d)
end

abstract type AbstractQGramMatchCounter end

abstract type AbstractQGramDistance <: SemiMetric end

function (dist::AbstractQGramDistance)(s1, s2)
	(s1 === missing) | (s2 === missing) && return missing
	counter = eval_start(dist)
	for (n1, n2) in _count(qgrams(s1, dist.q), qgrams(s2, dist.q))
		counter = eval_op(dist, counter, n1, n2)
	end
	eval_reduce(dist, counter)
end


"""
	QGram(q::Int)

Creates a QGram distance.

The distance corresponds to

``||v(s1, q) - v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the number of times a q-gram appears for the string s
"""
struct QGram <: AbstractQGramDistance
	q::Int
end
eval_start(::QGram) = 0
@inline eval_op(::QGram, c, n1::Integer, n2::Integer) = c + abs(n1 - n2)
eval_reduce(::QGram, c) = c

"""
	Cosine(q::Int)

Creates a Cosine distance.

The distance corresponds to

`` 1 - v(s1, q).v(s2, q)  / ||v(s1, q)|| * ||v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the  number of times a q-gram appears for the string s
"""
struct Cosine <: AbstractQGramDistance
	q::Int
end
eval_start(::Cosine) = (0, 0, 0)
@inline eval_op(::Cosine, c, n1::Integer, n2::Integer) = (c[1] + n1^2, c[2] + n2^2, c[3] + n1 * n2)
eval_reduce(::Cosine, c) = 1 - c[3] / sqrt(c[1] * c[2])

"""
	Jaccard(q::Int)

Creates a Jaccard distance.

The distance corresponds to 

``1 - |Q(s1, q) ∩ Q(s2, q)| / |Q(s1, q) ∪ Q(s2, q))|``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Jaccard <: AbstractQGramDistance
	q::Int
end
eval_start(::Jaccard) = (0, 0, 0)
@inline eval_op(::Jaccard, c, n1::Integer, n2::Integer) = (c[1] + (n1 > 0), c[2] + (n2 > 0), c[3] + (n1 > 0) * (n2 > 0))
eval_reduce(::Jaccard, c) = 1 - c[3] / (c[1] + c[2] - c[3])

"""
	SorensenDice(q::Int)

Creates a SorensenDice distance.

The distance corresponds to  

``1 - 2 * |Q(s1, q) ∩ Q(s2, q)|  / (|Q(s1, q)| + |Q(s2, q))|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct SorensenDice <: AbstractQGramDistance
	q::Int
end
eval_start(::SorensenDice) = (0, 0, 0)
@inline eval_op(::SorensenDice, c, n1::Integer, n2::Integer) = (c[1] + (n1 > 0), c[2] + (n2 > 0), c[3] + (n1 > 0) * (n2 > 0))
eval_reduce(::SorensenDice, c) = 1 - 2 * c[3] / (c[1] + c[2])

"""
	Overlap(q::Int)

Creates a Overlap distance.

The distance corresponds to  

``1 - |Q(s1, q) ∩ Q(s2, q)|  / min(|Q(s1, q)|, |Q(s2, q)|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Overlap <: AbstractQGramDistance
	q::Int
end
eval_start(::Overlap) = (0, 0, 0)
@inline eval_op(::Overlap, c, n1::Integer, n2::Integer) = (c[1] + (n1 > 0), c[2] + (n2 > 0), c[3] + (n1 > 0) * (n2 > 0))
eval_reduce(::Overlap, c) = 1 - c[3] / min(c[1], c[2])

"""
	NMD(q::Int)
	NMD(q::Int)

Creates a NMD (Normalized Multiset Distance) as introduced by Besiris and
Zigouris 2013. The goal with this distance is to behave similarly to a normalized
compression distance without having to do any actual compression (and thus being
faster to compute).

The distance corresponds to

``(sum(max.(m(s1), m(s2)) - min(M(s1), M(s2))) / max(M(s1), M(s2))``

where ``m(s)`` is the vector of q-gram counts for string ``s`` and ``M(s)`` is the
sum of those counts.

For details see:
https://www.sciencedirect.com/science/article/pii/S1047320313001417
"""
struct NMD <: AbstractQGramDistance
	q::Int
end
eval_start(::NMD) = (0, 0, 0)
@inline eval_op(::NMD, c, n1::Integer, n2::Integer) = (c[1] + n1, c[2] + n2, c[3] + max(n1, n2))
eval_reduce(::NMD, c) = (c[3] - min(c[1], c[2])) / max(c[1], c[2])

"""
	MorisitaOverlap(q::Int)

Creates a MorisitaOverlap distance, a general, statistical measure of
dispersion which can also be used on dictionaries such as created
from q-grams. See https://en.wikipedia.org/wiki/Morisita%27s_overlap_index
This is more fine-grained than many of the other QGramDistances since
it is based on the counts per q-gram rather than only which q-grams are
in the strings.

The distance corresponds to

``(2 * sum(m(s1) .* m(s2)) / (sum(m(s1).^2)*M(s2)/M(s1) + sum(m(s2).^2)*M(s1)/M(s2))``

where ``m(s)`` is the vector of q-gram counts for string ``s`` and ``M(s)`` is the
sum of those counts.
"""
struct MorisitaOverlap <: AbstractQGramDistance
	q::Int
end
eval_start(::MorisitaOverlap) = (0, 0, 0, 0, 0)
@inline eval_op(::MorisitaOverlap, c, n1::Integer, n2::Integer) = (c[1] + n1, c[2] + n2, c[3] + n1^2, c[4] + n2^2, c[5] + n1 * n2)
eval_reduce(::MorisitaOverlap, c) = 1 - 2 * c[5] / (c[3] * c[2] / c[1] + c[4] * c[1] / c[2])

