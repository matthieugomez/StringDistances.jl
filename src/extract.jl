"""
    extract(s1::AbstractString, iter, dist::PreMetric)

extrat returns the best element `iter` that has the best similarity score with `s1` according to the distance `dist`. 
The function is particularly fast for `Levenshtein` and `DamerauLevenshtein` distances (potentially modified by `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)
"""
function extract(s1::AbstractString, iter_s2, dist::Union{T, Partial{T}, TokenSort{T}, TokenSet{T}, TokenMax{T}}) where T <: Union{Levenshtein, DamerauLevenshtein}
    best_score = 0.0
    best_s2 = nothing
    for s2 in iter_s2
        score = compare(s1, s2, dist; min_dist = best_score)
        if (score !== missing) && (score > best_score)
            best_s2 = s2
            best_score = score
        end
    end
    return best_s2
end
function extract(s1::AbstractString, iter_s2, dist::PreMetric)
    best_score = 0.0
    best_s2 = nothing
    for s2 in iter_s2
        score = compare(s1, s2, dist)
        if (score !== missing) && (score > best_score)
            best_s2 = s2
            best_score = score
        end
    end
    return best_s2
end


function extract(::Missing, iter_s2, dist::PreMetric)
    return missing
end