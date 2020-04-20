"""
    findmax(s, itr, dist::StringDistance; min_score = 0.0) -> (x, index)

`findmax` returns the value and index of the element of `itr` that has the 
highest similarity score with `s` according to the distance `dist`. 
It returns `(nothing, nothing)` if none of the elements has a similarity score 
higher or equal to `min_score` (default to 0.0).

It is particularly optimized for [`Levenshtein`](@ref) and [`DamerauLevenshtein`](@ref) distances 
(as well as their modifications via [`Partial`](@ref), [`TokenSort`](@ref), [`TokenSet`](@ref), or [`TokenMax`](@ref)).

### Examples
```julia-repl
julia> using StringDistances
julia> s = "Newark"
julia> iter = ["New York", "Princeton", "San Francisco"]
julia> findmax(s, iter, Levenshtein())
("NewYork", 1)
julia> findmax(s, iter, Levenshtein(); min_score = 0.9)
(nothing, nothing)
```
"""
function Base.findmax(s, itr, dist::StringDistance; min_score = 0.0)
    min_score_atomic = Threads.Atomic{typeof(min_score)}(min_score)
    scores = [0.0 for _ in 1:Threads.nthreads()]
    is = [0 for _ in 1:Threads.nthreads()]
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, itr[i], dist; min_score = min_score_atomic[])
        score_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= score_old
            scores[Threads.threadid()] = score
            is[Threads.threadid()] = i
        end
    end
    imax = is[argmax(scores)]
    imax == 0 ? (nothing, nothing) : (itr[imax], imax)
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
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, itr[i], dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end
