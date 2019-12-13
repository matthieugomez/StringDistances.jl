"""
    findmax(s::AbstractString, itr, dist::StringDistance; min_score = 0.0) -> (x, index)

`findmax` returns the value and index of the element of `itr` that has the 
highest similarity score with `s` according to the distance `dist`. 
It returns `(nothing, nothing)` if none of the elements has a similarity score 
higher or equal to `min_score` (default to 0.0).
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances 
(as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`).
"""
function Base.findmax(s::AbstractString, itr, dist::StringDistance; min_score = 0.0)
    min_score = Threads.Atomic{typeof(min_score)}(min_score)
    scores = [0.0 for _ in 1:Threads.nthreads()]
    is = [0 for _ in 1:Threads.nthreads()]
    Threads.@threads for i in collect(keys(itr))
        score = compare(s, itr[i], dist; min_score = min_score[])
        score_old = Threads.atomic_max!(min_score, score)
        if score >= score_old
            scores[Threads.threadid()] = score
            is[Threads.threadid()] = i
        end
    end
    imax = is[argmax(scores)]
    imax == 0 ? (nothing, nothing) : (itr[imax], imax)
end

"""
    findall(s::AbstractString, itr, dist::StringDistance; min_score = 0.8)
    
`findall` returns the vector of indices for elements of `itr` that have a 
similarity score higher or equal than `min_score` according to the distance `dist`.
If there are no such elements, return an empty array. 
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances 
(as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`).
"""
function Base.findall(s::AbstractString, itr, dist::StringDistance; min_score = 0.8)
    out = [Int[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in collect(keys(itr))
        score = compare(s, itr[i], dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end
