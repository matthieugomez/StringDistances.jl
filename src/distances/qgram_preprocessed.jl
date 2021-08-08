# sometimes we already preprocess the strings
# We now define special methods for these special string types
"""
	QGramDict(s, q::Integer = 2)

An iterator with a pre-computed dictionary of its qgrams. This enables faster calculation of QGram 
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
struct QGramDict{S, K}
    s::S
    q::Int
    counts::Dict{K,Int}
end
Base.length(s::QGramDict) = length(s.s)
Base.iterate(s::QGramDict) = iterate(s.s)
Base.iterate(s::QGramDict, state) = iterate(s.s, state)

function QGramDict(s, q::Integer = 2)
    (s isa QGramDict) && (s.q == q) && return s
    qgs = qgrams(s, q)
    QGramDict{typeof(s), eltype(qgs)}(s, q, countdict(qgs))
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

function (dist::AbstractQGramDistance)(qc1::QGramDict, qc2::QGramDict)
    dist.q == qc1.q == qc2.q || throw(ArgumentError("The distance and the QGramDict must have the same qgram length"))
    counter = newcounter(dist)
    d1, d2 = qc1.counts, qc2.counts
    for (k1, c1) in d1
        index = Base.ht_keyindex2!(d2, k1)
		if index > 0
			count!(dist, counter, c1, d2.vals[index])
		else
			count!(dist, counter, c1, 0)
        end
    end
    for (k2, c2) in d2
        index = Base.ht_keyindex2!(d1, k2)
		if index <= 0
			count!(dist, counter, 0, c2)
        end
    end
    calculate(dist, counter)
end

"""
	QGramSortedVector(s, q::Integer = 2)

An iterator with a pre-computed sorted vector of its qgrams. This enables faster calculation of QGram 
distances.

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
struct QGramSortedVector{S, K}
    s::S
    q::Int
    counts::Vector{Pair{K,Int}}
end
Base.length(s::QGramSortedVector) = length(s.s)
Base.iterate(s::QGramSortedVector) = iterate(s.s)
Base.iterate(s::QGramSortedVector, state) = iterate(s.s, state)

function QGramSortedVector(s, q::Integer = 2)
    (s isa QGramSortedVector) && (s.q == q) && return s
    qgs = qgrams(s, q)
    countpairs = collect(countdict(qgs))
    sort!(countpairs, by = first)
    QGramSortedVector{typeof(s), eltype(qgs)}(s, q, countpairs)
end



# To implement the distances we will count qgram matches
# between strings or pre-calculated AbstractQgramCounts objects.
# The abstract type defines different fallback versions which can be
# specialied by subtypes for best performance.
function (dist::AbstractQGramDistance)(qc1::QGramSortedVector, qc2::QGramSortedVector)
    dist.q == qc1.q == qc2.q || throw(ArgumentError("The distance and the QGramSortedVectors must have the same qgram length"))
    counter = newcounter(dist)
    d1, d2 = qc1.counts, qc2.counts
    i1 = i2 = 1
    while true
    	# length can be zero
        if i2 > length(d2)
			for i in i1:length(d1)
				@inbounds count!(dist, counter, d1[i][2], 0)
            end
            break
        elseif i1 > length(d1)
			for i in i2:length(d2)
				@inbounds count!(dist, counter, 0, d2[i][2])
            end
            break
        end
        @inbounds k1, n1 = d1[i1]
        @inbounds k2, n2 = d2[i2]
        cmpval = Base.cmp(k1, k2)
		if cmpval == -1 # k1 < k2
			count!(dist, counter, n1, 0)
            i1 += 1
        elseif cmpval == +1 # k2 < k1
        	count!(dist, counter, 0, n2)
            i2 += 1
		else
			count!(dist, counter, n1, n2)
            i1 += 1
            i2 += 1
        end
    end
    calculate(dist, counter)
end

