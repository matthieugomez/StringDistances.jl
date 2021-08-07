# sometimes we already preprocess the strings as AbstractQGramCounts
# We know define how QgramDistances can be computed from these AbstractQGramCounts


abstract type AbstractQGramCounts{Q,K} end
q(qc::AbstractQGramCounts{Q,K}) where {Q,K} = Q
counts(qc::AbstractQGramCounts) = qc.counts
Base.length(qc::AbstractQGramCounts{Q}) where Q = length(qc.counts) + Q - 1

function (dist::AbstractQGramDistance)(qc1::QC, qc2::QC) where {QC<:AbstractQGramCounts}
    @assert dist.q == q(qc1)
	@assert dist.q == q(qc2)
	counter = newcounter(dist)
	countmatches!(counter, counts(qc1), counts(qc2))
    calculate(dist, counter)
end

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

QGramDict(s, q::Integer = 2) = QGramDict(collect(s), q)


function countmatches!(mc::AbstractQGramMatchCounter, d1::Dict{K,I}, d2::Dict{K,I}) where {K,I<:Integer}
    for (k1, c1) in d1
        index = Base.ht_keyindex2!(d2, k1)
		if index > 0
			count!(mc, c1, d2.vals[index])
		else
			count!(mc, c1, 0)
        end
    end
    for (k2, c2) in d2
        index = Base.ht_keyindex2!(d1, k2)
		if index <= 0
			count!(mc, 0, c2)
        end
    end
end


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
function countmatches!(mc::AbstractQGramMatchCounter, d1::Vector{Pair{K,I}}, d2::Vector{Pair{K,I}}) where {K,I<:Integer}
    i1 = i2 = 1
    while true
    	# length can be zero
        if i2 > length(d2)
			for i in i1:length(d1)
				@inbounds count!(mc, d1[i][2], 0)
            end
            return
        elseif i1 > length(d1)
			for i in i2:length(d2)
				@inbounds count!(mc, 0, d2[i][2])
            end
            return
        end
        @inbounds k1, n1 = d1[i1]
        @inbounds k2, n2 = d2[i2]
        cmpval = Base.cmp(k1, k2)
		if cmpval == -1 # k1 < k2
			count!(mc, n1, 0)
            i1 += 1
        elseif cmpval == +1 # k2 < k1
        	count!(mc, 0, n2)
            i2 += 1
		else
			count!(mc, n1, n2)
            i1 += 1
            i2 += 1
        end
    end
end


