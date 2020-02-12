struct QGramIterator{S <: Union{AbstractString, AbstractVector}}
	s::S   # Collection
	q::Int # Length of Qgram
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


#q-grams of AbstractVector
# Alternatively, I could also use partition in IterTools but it creates a vector for each iteration
# so it does not seem to be worth it.
function Base.iterate(qgram::QGramIterator{<: AbstractVector}, state = firstindex(qgram.s))
	state + qgram.q - 1 > lastindex(qgram.s) && return nothing
	view(qgram.s, state:(state + qgram.q - 1)), state + 1
end
Base.eltype(qgram::QGramIterator{<: AbstractVector}) = typeof(first(qgram))


"""
Return an iterator on the q-gram of a string

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
qgrams(s::Union{AbstractString, AbstractVector}, q::Integer) = QGramIterator(s, q)
qgrams(s, q::Integer) = QGramIterator(collect(s), q)


# For two iterators x1 and x2, that define a length and eltype method,
# this returns a dictionary which, for each element in x1 or x2, 
# returns a tuple with the numbers of times it appears in x1 and x2
function count_map(s1, s2)
	K = promote_type(eltype(s1), eltype(s2))
	d = Dict{K, Tuple{Int, Int}}()
	sizehint!(d, length(s1) + length(s2))
	# I use a faster way to change a dictionary key
	# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
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


abstract type QGramDistance <: SemiMetric end

"""
	QGram(q::Int)

Creates a QGram metric.

The distance corresponds to

``||v(s1, q) - v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the number of times a q-gram appears for the string s
"""
struct QGram <: QGramDistance
	q::Int
end

function (dist::QGram)(s1, s2)
	(ismissing(s1) | ismissing(s2)) && return missing
	itr = values(count_map(qgrams(s1, dist.q), qgrams(s2, dist.q)))
	n = 0
	for (n1, n2) in itr
		n += abs(n1 - n2)
	end
	n
end

"""
	Cosine(q::Int)

Creates a Cosine metric.

The distance corresponds to

`` 1 - v(s1, q).v(s2, q)  / ||v(s1, q)|| * ||v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the  number of times a q-gram appears for the string s
"""
struct Cosine <: QGramDistance
	q::Int
end

function (dist::Cosine)(s1, s2, max_dist = nothing)
	(ismissing(s1) | ismissing(s2)) && return missing
	itr = values(count_map(qgrams(s1, dist.q), qgrams(s2, dist.q)))
	norm1, norm2, prodnorm = 0, 0, 0
	for (n1, n2) in itr
		norm1 += n1^2
		norm2 += n2^2
		prodnorm += n1 * n2
	end
	1.0 - prodnorm / (sqrt(norm1) * sqrt(norm2))
end

"""
	Jaccard(q::Int)

Creates a Jaccard metric.

The distance corresponds to 

``1 - |Q(s1, q) ∩ Q(s2, q)| / |Q(s1, q) ∪ Q(s2, q))|``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Jaccard <: QGramDistance
	q::Int
end

function (dist::Jaccard)(s1, s2, max_dist = nothing)
	(ismissing(s1) | ismissing(s2)) && return missing
	itr = values(count_map(qgrams(s1, dist.q), qgrams(s2, dist.q)))
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in itr
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / (ndistinct1 + ndistinct2 - nintersect)
end

"""
	SorensenDice(q::Int)

Creates a SorensenDice metric

The distance corresponds to  

``1 - 2 * |Q(s1, q) ∩ Q(s2, q)|  / (|Q(s1, q)| + |Q(s2, q))|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct SorensenDice <: QGramDistance
	q::Int
end

function (dist::SorensenDice)(s1, s2, max_dist = nothing)
	(ismissing(s1) | ismissing(s2)) && return missing
	itr = values(count_map(qgrams(s1, dist.q), qgrams(s2, dist.q)))
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in itr
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - 2.0 * nintersect / (ndistinct1 + ndistinct2)
end

"""
	Overlap(q::Int)

Creates a Overlap metric

The distance corresponds to  

``1 - |Q(s1, q) ∩ Q(s2, q)|  / min(|Q(s1, q)|, |Q(s2, q)|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Overlap <: QGramDistance
	q::Int
end

function (dist::Overlap)(s1, s2, max_dist = nothing)
	(ismissing(s1) | ismissing(s2)) && return missing
	itr = values(count_map(qgrams(s1, dist.q), qgrams(s2, dist.q)))
	ndistinct1, ndistinct2, nintersect = 0, 0, 0
	for (n1, n2) in itr
		ndistinct1 += n1 > 0
		ndistinct2 += n2 > 0
		nintersect += (n1 > 0) & (n2 > 0)
	end
	1.0 - nintersect / min(ndistinct1, ndistinct2)
end
