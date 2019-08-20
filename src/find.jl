"""
    find_best(s1::AbstractString, iter, dist::PreMetric)

`find_best` returns the element of the iterator `iter` that has the highest similarity score with `s1` according to the distance `dist`. 
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function find_best(s1::AbstractString, iter_s2, dist::Union{T, Partial{T}, TokenSort{T}, TokenSet{T}, TokenMax{T}}) where T <: Union{Levenshtein, DamerauLevenshtein}
    best_score = 0.0
    best_s2 = nothing
    for s2 in iter_s2
        score = compare(s1, s2, dist; min_score = best_score)
        if score > best_score
            best_s2 = s2
            best_score = score
        end
    end
    return best_s2
end
function find_best(s1::AbstractString, iter_s2, dist::PreMetric)
    best_score = 0.0
    best_s2 = nothing
    for s2 in iter_s2
        score = compare(s1, s2, dist)
        if score > best_score
            best_s2 = s2
            best_score = score
        end
    end
    return best_s2
end


"""
    find_all(s1::AbstractString, iter, dist::PreMetric; min_score = 0.8)
`find_all` returns the vector with all the elements of `iter` that have a similarity score higher than `min_score` according to the distance `dist`. 
The function is optimized for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function find_all(s1::AbstractString, iter_s2, dist::Union{T, Partial{T}, TokenSort{T}, TokenSet{T}, TokenMax{T}}; min_score = 0.8) where T <: Union{Levenshtein, DamerauLevenshtein}
    collect(s2 for s2 in iter_s2 if compare(s1, s2, dist; min_score = min_score) >= min_score)
end
function find_all(s1::AbstractString, iter_s2, dist::PreMetric; min_score = 0.8)
    collect(s2 for s2 in iter_s2 if compare(s1, s2, dist) >= min_score)
end
