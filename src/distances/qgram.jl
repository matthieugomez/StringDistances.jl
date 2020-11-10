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
function _count(s1, s2)
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
	return values(d)
end

# Turn a sequence of qgrams to a count dict for them, i.e. map each
# qgram to the number of times it has been seen.
function countdict(qgrams)
    d = Dict{eltype(qgrams), Int32}()
    for qg in qgrams
        index = Base.ht_keyindex2!(d, qg)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = qg
			@inbounds d.vals[index] = d.vals[index][1] + 1
		else
			@inbounds Base._setindex!(d, 1, qg, -index)
		end
    end
    d
end

abstract type AbstractQGramCounts{Q,K} end
q(qc::AbstractQGramCounts{Q,K}) where {Q,K} = Q
counts(qc::AbstractQGramCounts) = qc.counts
Base.length(qc::AbstractQGramCounts{Q}) where Q = length(qc.counts) + Q - 1
"""
	QGramDict(s, q::Integer = 2)

Creates a QGramDict that pre-calculates (pre-counts) the qgrams
of a string or stream. This enables faster calculation of QGram 
distances.

Note that the qgram length must correspond with the q length used
in the distance.

## Examples
```julia
str1, str2 = "my string", "another string"
qd1 = QGramDict(str1, 2)
qd2 = QGramDict(str2, 2)
evaluate(Overlap(2), qd1, qd2)
```
"""
struct QGramDict{Q,K} <: AbstractQGramCounts{Q,K}
    counts::Dict{K,Int}
end
function QGramDict(s::Union{AbstractString, AbstractVector}, q::Integer = 2)
    @assert q >= 1
    qgs = qgrams(s, q)
    QGramDict{q, eltype(qgs)}(countdict(qgs))
end
QGramDict(s, q::Integer = 2) = QGramDict(collect(s), q)

"""
	QGramSortedVector(s, q::Integer = 2)

Creates a QGramSortedVector that pre-calculates (pre-counts) the 
qgrams of a string or stream. This enables faster calculation of
QGram distances.

Since qgrams are sorted in lexicographic order QGram distances can be 
calculated even faster than when using a QGramDict. However, the 
sorting means that updating the counts after creation is less 
efficient. However, for most use cases QGramSortedVector is preferred
over a QgramDict.

Note that the qgram length must correspond with the q length used
in the distance.

## Examples
```julia
str1, str2 = "my string", "another string"
qs1 = QGramSortedVector(str1, 2)
qs2 = QGramSortedVector(str2, 2)
evaluate(Jaccard(2), qs1, qs2)
```
"""
struct QGramSortedVector{Q,K} <: AbstractQGramCounts{Q,K}
    counts::Vector{Pair{K,Int}}
end
function QGramSortedVector(s::Union{AbstractString, AbstractVector}, q::Integer = 2)
    @assert q >= 1
    qgs = qgrams(s, q)
    countpairs = collect(countdict(qgs))
    sort!(countpairs, by = first)
    QGramSortedVector{q, eltype(qgs)}(countpairs)
end
QGramSortedVector(s, q::Integer = 2) = QGramSortedVector(collect(s), q)

# To implement the distances we will count qgram matches
# between strings or pre-calculated AbstractQgramCounts objects.
# The abstract type defines different fallback versions which can be
# specialied by subtypes for best performance.
abstract type AbstractQGramMatchCounter end
@inline countleft!(c::AbstractQGramMatchCounter, qg, n1::Integer) = countleft!(c, n1)
@inline countright!(c::AbstractQGramMatchCounter, qg, n2::Integer) = countright!(c, n2)
@inline countboth!(c::AbstractQGramMatchCounter, qg, n1::Integer, n2::Integer) =
	countboth!(c, n1, n2)
@inline function countboth!(c::AbstractQGramMatchCounter, n1::Integer, n2::Integer)
	countleft!(c, n1)
	countright!(c, n2)
	countshared!(c, n1, n2)
end
@inline countshared!(c::AbstractQGramMatchCounter, qg, n1::Integer, n2::Integer) = countshared!(c, n1, n2)

# Subtypes must implement these methods:
@inline countleft!(c::AbstractQGramMatchCounter, n1::Integer) =
	error("countleft! not implemented for $(typeof(c))")
@inline countright!(c::AbstractQGramMatchCounter, n2::Integer) =
	error("countright! not implemented for $(typeof(c))")

# Subtypes either must overwrite countboth! from above (so it not uses countshared!) or implement:
@inline countshared!(c::AbstractQGramMatchCounter, n1::Integer, n2::Integer) =
	error("countshared! not implemented for $(typeof(c))")

function countmatches!(mc::AbstractQGramMatchCounter, d1::Vector{Pair{K,I}}, d2::Vector{Pair{K,I}}) where {K,I<:Integer}
    i1 = i2 = 1
    while i1 <= length(d1) || i2 <= length(d2)
        if i2 > length(d2)
			for i in i1:length(d1)
				@inbounds countleft!(mc, d1[i][1], d1[i][2])
            end
            return
        elseif i1 > length(d1)
			for i in i2:length(d2)
				@inbounds countright!(mc, d2[i][1], d2[i][2])
            end
            return
        end
        @inbounds k1, n1 = d1[i1]
        @inbounds k2, n2 = d2[i2]
        cmpval = Base.cmp(k1, k2)
		if cmpval == -1 # k1 < k2
			countleft!(mc, k1, n1)
            i1 += 1
        elseif cmpval == +1 # k2 < k1
			countright!(mc, k2, n2)
            i2 += 1
		else
			countboth!(mc, k1, n1, n2)
            i1 += 1
            i2 += 1
        end
    end
end

function countmatches!(mc::AbstractQGramMatchCounter, d1::Dict{K,I}, d2::Dict{K,I}) where {K,I<:Integer}
    for (k1, c1) in d1
        index = Base.ht_keyindex2!(d2, k1)
		if index > 0
			countboth!(mc, k1, c1, d2.vals[index])
		else
			countleft!(mc, k1, c1)
        end
    end
    for (k2, c2) in d2
        index = Base.ht_keyindex2!(d1, k2)
		if index <= 0
			countright!(mc, k2, c2)
        end
    end
end

abstract type QGramDistance <: SemiMetric end

function (dist::QGramDistance)(s1, s2)
	((s1 === missing) | (s2 === missing)) && return missing
	counter = newcounter(dist)
	for (n1, n2) in _count(qgrams(s1, dist.q), qgrams(s2, dist.q))
		countboth!(counter, n1, n2)
	end
	calculate(dist, counter)
end

function (dist::QGramDistance)(qc1::QC, qc2::QC) where {QC<:AbstractQGramCounts}
    @assert dist.q == q(qc1)
	@assert dist.q == q(qc2)
	counter = newcounter(dist)
	countmatches!(counter, counts(qc1), counts(qc2))
    calculate(dist, counter)
end

"""
	QGram(q::Int)

Creates a QGram distance.

The distance corresponds to

``||v(s1, q) - v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the number of times a q-gram appears for the string s
"""
struct QGram <: QGramDistance
	q::Int
end

mutable struct SingleCounter{T, QD<:QGramDistance} <: AbstractQGramMatchCounter
	n::T
end

newcounter(d::QGram) = SingleCounter{Int, QGram}(0)

@inline countleft!(c::SingleCounter{Int, QGram}, n1::Integer) = c.n += n1 # n1 === abs(n1 - 0)
@inline countright!(c::SingleCounter{Int, QGram}, n2::Integer) = c.n += n2 # n2 === abs(0 - n2)
@inline countboth!(c::SingleCounter{Int, QGram}, n1::Integer, n2::Integer) = c.n += abs(n1 - n2)

calculate(dist::QGram, c::SingleCounter{Int, QGram}) = c.n

"""
	Cosine(q::Int)

Creates a Cosine distance.

The distance corresponds to

`` 1 - v(s1, q).v(s2, q)  / ||v(s1, q)|| * ||v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the  number of times a q-gram appears for the string s
"""
struct Cosine <: QGramDistance
	q::Int
end

mutable struct ThreeCounters{T, QD<:QGramDistance} <: AbstractQGramMatchCounter
	left::T
	right::T
	shared::T
end

newcounter(d::Cosine) = ThreeCounters{Int, Cosine}(0, 0, 0)

@inline countleft!(c::ThreeCounters{Int, Cosine}, n1::Integer) = c.left += n1^2
@inline countright!(c::ThreeCounters{Int, Cosine}, n2::Integer) = c.right += n2^2
@inline countshared!(c::ThreeCounters{Int, Cosine}, n1::Integer, n2::Integer) = c.shared += n1 * n2

calculate(d::Cosine, c::ThreeCounters{Int, Cosine}) =
	1.0 - c.shared / (sqrt(c.left) * sqrt(c.right))

"""
	Jaccard(q::Int)

Creates a Jaccard distance.

The distance corresponds to 

``1 - |Q(s1, q) ∩ Q(s2, q)| / |Q(s1, q) ∪ Q(s2, q))|``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Jaccard <: QGramDistance
	q::Int
end

calculate(d::Jaccard, c::ThreeCounters{Int, Jaccard}) =
	1.0 - c.shared / (c.left + c.right - c.shared)

"""
	SorensenDice(q::Int)

Creates a SorensenDice distance.

The distance corresponds to  

``1 - 2 * |Q(s1, q) ∩ Q(s2, q)|  / (|Q(s1, q)| + |Q(s2, q))|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct SorensenDice <: QGramDistance
	q::Int
end

calculate(d::SorensenDice, c::ThreeCounters{Int, SorensenDice}) =
	1.0 - 2.0 * c.shared / (c.left + c.right)

"""
	Overlap(q::Int)

Creates a Overlap distance.

The distance corresponds to  

``1 - |Q(s1, q) ∩ Q(s2, q)|  / min(|Q(s1, q)|, |Q(s2, q)|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Overlap <: QGramDistance
	q::Int
end

const IntersectionDist = Union{Jaccard, SorensenDice, Overlap}
newcounter(d::IntersectionDist) = ThreeCounters{Int, typeof(d)}(0, 0, 0)

@inline countleft!(c::ThreeCounters{Int, QD}, n1::Integer) where {QD<:IntersectionDist} =
	c.left += (n1 > 0)
@inline countright!(c::ThreeCounters{Int, QD}, n2::Integer) where {QD<:IntersectionDist} =
	c.right += (n2 > 0)
@inline countshared!(c::ThreeCounters{Int, QD}, n1::Integer, n2::Integer) where {QD<:IntersectionDist} =
	c.shared += (n1 > 0) & (n2 > 0)

calculate(d::Overlap, c::ThreeCounters{Int, Overlap}) =
	1.0 - c.shared / min(c.left, c.right)

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
struct MorisitaOverlap <: QGramDistance
	q::Int
end

mutable struct FiveCounters{T, QD<:QGramDistance} <: AbstractQGramMatchCounter
	leftsum::T    # sum(m(s1))
	rightsum::T   # sum(m(s2))
	leftsq::T     # sum(m(s1).^2)
	rightsq::T    # sum(m(s2).^2)
	shared::T     # sum(m(s1) .* m(s2))
end

newcounter(d::MorisitaOverlap) = FiveCounters{Int, MorisitaOverlap}(0, 0, 0, 0, 0)

@inline function countleft!(c::FiveCounters{Int, MorisitaOverlap}, n1::Integer)
	c.leftsum += n1
	c.leftsq += (n1^2)
end

@inline function countright!(c::FiveCounters{Int, MorisitaOverlap}, n2::Integer)
	c.rightsum += n2
	c.rightsq += (n2^2)
end

@inline countshared!(c::FiveCounters{Int, MorisitaOverlap}, n1::Integer, n2::Integer) =
	c.shared += (n1 * n2)

calculate(d::MorisitaOverlap, c::FiveCounters{Int, MorisitaOverlap}) =
	1.0 - ((2 * c.shared) / (c.leftsq*c.rightsum/c.leftsum + c.rightsq*c.leftsum/c.rightsum))

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
struct NMD <: QGramDistance
	q::Int
end

newcounter(d::NMD) = ThreeCounters{Int, NMD}(0, 0, 0)

@inline function countleft!(c::ThreeCounters{Int, NMD}, n1::Integer)
	c.left += n1
	c.shared += n1 # max(n1, 0) == n1
end

@inline function countright!(c::ThreeCounters{Int, NMD}, n2::Integer)
	c.right += n2
	c.shared += n2 # max(n2, 0) == n2
end

@inline function countboth!(c::ThreeCounters{Int, NMD}, n1::Integer, n2::Integer)
	c.left += n1
	c.right += n2
	c.shared += max(n1, n2)
end

calculate(d::NMD, c::ThreeCounters{Int, NMD}) =
	(c.shared - min(c.left, c.right)) / max(c.left, c.right)
