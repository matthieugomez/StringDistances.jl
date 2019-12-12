"""
    findmax(s::AbstractString, iter::AbstractVector, dist::StringDistance; min_score = 0.0)

`findmax` returns the value and index of the element of `iter` that has the highest similarity score with `s` according to the distance `dist`. 
It returns `(nothing, nothing)` if none of the elements has a similarity score higher or equal to `min_score` (default to 0.0)
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function Base.findmax(s::AbstractString, iter::AbstractVector, dist::StringDistance; min_score = 0.0)
    min_score >= 0 || throw("min_score should be positive")
    is = [0 for _ in 1:Threads.nthreads()]
    xs = ["" for _ in 1:Threads.nthreads()]
    scores = [-1.0 for _ in 1:Threads.nthreads()]
    min_score_atomic = Threads.Atomic{typeof(min_score)}(min_score)
    Threads.@threads for i in 1:length(iter)
        score = compare(s, iter[i], dist; min_score = min_score_atomic[])
        min_score_atomic_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= min_score_atomic_old
            score == 1.0 && return i
            is[Threads.threadid()] = i
            xs[Threads.threadid()] = iter[i]
            scores[Threads.threadid()] = score
        end
    end
    i = argmax(scores)
    is[i] == 0 ? (nothing, nothing) : (xs[i], is[i])
end


"""
    findall(s::AbstractString, iter::AbstractVector, dist::StringDistance; min_score = 0.8)
`findall` returns the vector of indices for elements of `iter` that have a similarity score higher or equal than `min_score` according to the distance `dist`. 
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function Base.findall(s::AbstractString, iter::AbstractVector, dist::StringDistance; min_score = 0.8)
    out = [Int[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:length(iter)
        score = compare(s, iter[i], dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end
