"""
    findnearest(s, itr, dist::StringDistance) -> (x, index)

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
function findnearest(s, itr, dist::StringDistance; min_score = 0.0)
    min_score_atomic = Threads.Atomic{Float64}(min_score)
    scores = [0.0 for _ in 1:Threads.nthreads()]
    is = [0 for _ in 1:Threads.nthreads()]
    s = _preprocess(dist, s)
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, _preprocess(dist, itr[i]), dist; min_score = min_score_atomic[])
        score_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= score_old
            scores[Threads.threadid()] = score
            is[Threads.threadid()] = i
        end
    end
    imax = is[argmax(scores)]
    imax == 0 ? (nothing, nothing) : (itr[imax], imax)
end
_preprocess(dist::AbstractQGramDistance, ::Missing) = missing
_preprocess(dist::AbstractQGramDistance, s) = QGramSortedVector(s, dist.q)
_preprocess(dist::StringDistance, s) = s

function Base.findmax(s, itr, dist::StringDistance; min_score = 0.0)
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
function Base.findall(s, itr, dist::StringDistance; min_score = 0.8)
    out = [Int[] for _ in 1:Threads.nthreads()]
    s = _preprocess(dist, s)
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, _preprocess(dist, itr[i]), dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end
