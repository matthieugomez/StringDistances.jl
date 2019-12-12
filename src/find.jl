"""
    find_best(s1::AbstractString, iter, dist::PreMetric; min_score = 0.0)

`find_best` returns the element of the iterator `iter` that has the highest similarity score with `s1` according to the distance `dist`. Return nothing if all elements have a similarity score below `min_score`.
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function find_best(s1::AbstractString, iter_s2, dist::PreMetric; min_score = 0.0)
    min_score >= 0 || throw("min_score should be positive")
    best_s2s = AbstractString["" for _ in 1:Threads.nthreads()]
    best_scores = [-1.0 for _ in 1:Threads.nthreads()]
    min_score_atomic = Threads.Atomic{typeof(min_score)}(min_score)
    Threads.@threads for s2 in iter_s2
        score = compare(s1, s2, dist; min_score = min_score_atomic[])
        min_score_atomic_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= min_score_atomic_old
            best_s2s[Threads.threadid()] = s2
            best_scores[Threads.threadid()] = score
            score == 1.0 && return s2
        end
    end
    i = argmax(best_scores)
    if best_scores[i] < 0
        return nothing
    else
        return best_s2s[i]
    end
end


"""
    find_all(s1::AbstractString, iter, dist::PreMetric; min_score = 0.8)
`find_all` returns the vector with all the elements of `iter` that have a similarity score higher or equal than `min_score` according to the distance `dist`. 
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function find_all(s1::AbstractString, iter_s2, dist::PreMetric; min_score = 0.8)
    best_s2s = [eltype(iter_s2)[] for _ in 1:Threads.nthreads()]
    Threads.@threads for s2 in iter_s2
        score = compare(s1, s2, dist; min_score = min_score)
        if score >= min_score
            push!(best_s2s[Threads.threadid()], s2)
        end
    end
    vcat(best_s2s...)
end
