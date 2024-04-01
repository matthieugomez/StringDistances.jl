"""
    compare(s1, s2, dist)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the distance `dist`.

### Examples
```julia-repl
julia> compare("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function compare(s1, s2, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    1 - Normalized(dist)(s1, s2; max_dist = 1 - min_score)
end 

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
julia> iter = ["New York", "Princeton", "San Francisco"]
julia> findnearest(s, iter, Levenshtein())
("NewYork", 1)
julia> findnearest(s, iter, Levenshtein(); min_score = 0.9)
(nothing, nothing)
```
"""
function findnearest(s, itr, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    scores = tmap(itr) do x
        compare(s, _preprocess(dist, x), dist; min_score = min_score)
    end

    imax = argmax(scores)
    iszero(scores) ? (nothing, nothing) : (itr[imax], imax)
end

_preprocess(dist::AbstractQGramDistance, ::Missing) = missing
_preprocess(dist::AbstractQGramDistance, s) = QGramSortedVector(s, dist.q)
_preprocess(dist::Union{StringSemiMetric, StringMetric}, s) = s

function Base.findmax(s, itr, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    @warn "findmax(s, itr, dist; min_score) is deprecated. Use findnearest(s, itr, dist; min_score)"
    findnearest(s, itr, dist; min_score = min_score)
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
    out = Int[]
    s = _preprocess(dist, s)
    @tasks for i in eachindex(itr)
        score = compare(s, _preprocess(dist, itr[i]), dist; min_score = min_score)
        @one_by_one if score >= min_score
            push!(out, i)
        end
    end
    return out
end
