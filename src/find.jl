"""
    similarity(s1, s2, dist)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the distance `dist`.

### Examples
```julia-repl
julia> similarity("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function similarity(s1, s2, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    1 - Normalized(dist)(s1, s2; max_dist = 1 - min_score)
end

Base.@deprecate compare(s1, s2, dist; min_score = 0.0) similarity(s1, s2, dist; min_score = min_score)

"""
    findnearest(s, itr, dist::Union{StringMetric, StringSemiMetric}) -> (x, index)

`findnearest` returns the value and index of the element of `itr` that has the 
lowest distance with `s` according to the distance `dist`. 

It is particularly optimized for [`Levenshtein`](@ref) and [`DamerauLevenshtein`](@ref) distances 
(as well as their modifications via [`Partial`](@ref), [`TokenSort`](@ref), [`TokenSet`](@ref), or [`TokenMax`](@ref)).

### Examples
```julia-repl
julia> using StringDistances
julia> s = "Newark"
julia> iter = ["NewYork", "Princeton", "San Francisco"]
julia> findnearest(s, iter, Levenshtein())
("NewYork", 1)
julia> findnearest(s, iter, Levenshtein(); min_score = 0.9)
(nothing, nothing)
```
"""
function findnearest(s, itr, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    _citr = collect(itr)
    isempty(_citr) && return (nothing, nothing)

    _preprocessed_s = _preprocess(dist, s)
    _preprocessed_s === missing && return (nothing, nothing)
    ranges = _chunk_ranges(length(_citr))
    scores = Vector{Float64}(undef, length(_citr))

    Threads.@threads :dynamic for ir in eachindex(ranges)
        for i in ranges[ir]
            scores[i] = _similarity_or_neg_inf(_preprocessed_s, _citr[i], dist; min_score = min_score)
        end
    end

    imax = argmax(scores)
    scores[imax] > 0 ? (_citr[imax], imax) : (nothing, nothing)
end

_preprocess(dist::AbstractQGramDistance, ::Missing) = missing
_preprocess(dist::AbstractQGramDistance, s) = QGramSortedVector(s, dist.q)
_preprocess(dist::Union{StringSemiMetric, StringMetric}, s) = s

function _similarity_or_neg_inf(_preprocessed_s, x, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    x === missing && return -Inf
    _preprocessed_x = _preprocess(dist, x)
    _preprocessed_x === missing && return -Inf
    score = similarity(_preprocessed_s, _preprocessed_x, dist; min_score = min_score)
    ismissing(score) ? -Inf : score
end


"""
    findall(s, itr , dist::StringDistance; min_score = 0.8)
    
`findall` returns the vector of indices for elements of `itr` that have a 
similarity score higher or equal than `min_score` according to the distance `dist`.
If there are no such elements, return an empty array. 

It is particularly optimized for [`Levenshtein`](@ref) and [`DamerauLevenshtein`](@ref) distances 
(as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`).

### Examples
```julia-repl
julia> using StringDistances
julia> s = "Newark"
julia> iter = ["Newwark", "Princeton", "San Francisco"]
julia> findall(s, iter, Levenshtein())
1-element Array{Int64,1}:
 1
julia> findall(s, iter, Levenshtein(); min_score = 0.9)
0-element Array{Int64,1}
```
"""
function Base.findall(s, itr, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.8)
    _citr = collect(itr)
    isempty(_citr) && return Int[]
    _preprocessed_s = _preprocess(dist, s)
    _preprocessed_s === missing && return Int[]

    ranges = _chunk_ranges(length(_citr))
    scores = Vector{Float64}(undef, length(_citr))

    Threads.@threads :dynamic for ir in eachindex(ranges)
        for i in ranges[ir]
            scores[i] = _similarity_or_neg_inf(_preprocessed_s, _citr[i], dist; min_score = min_score)
        end
    end

    return findall(>=(min_score), scores)
end
