
struct Normalized{S <: SemiMetric} <: SemiMetric
    dist::S
end
# A normalized distance is between 0 and 1, and accept a third argument, max_dist.
function (dist::Normalized{<: Union{Levenshtein, DamerauLevenshtein}})(s1, s2, max_dist = 1.0)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    d = dist.dist(s1, s2, ceil(Int, len2 * max_dist))
    out = d / len2
    out > max_dist ? 1.0 : out
end

function (dist::Normalized{<: QGramDistance})(s1, s2, max_dist = 1.0)
    ((s1 === missing) | (s2 === missing)) && return missing
    # When string length < q for qgram distance, returns s1 == s2
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 <= dist.dist.q - 1 && return convert(Float64, s1 != s2)
    if typeof(dist.dist) <: QGram
        dist.dist(s1, s2) / (len1 + len2 - 2 * dist.dist.q + 2)
    else
        dist.dist(s1, s2)
    end
end

function (dist::Normalized)(s1, s2, max_dist = 1.0)
    dist.dist(s1, s2)
end

"""
   normalize(dist::SemiMetric)

   Normalize a metric, so that `evaluate` always return a Float64 between 0 and 1
"""
normalize(dist::SemiMetric) = Normalized(dist)
normalize(dist::Normalized) = dist



"""
   Winkler(dist; p::Real = 0.1, threshold::Real = 0.7, maxlength::Integer = 4)

Creates the `Winkler{dist, p, threshold, maxlength}` distance.

`Winkler{dist, p, threshold, length)` modifies the string distance `normalize(dist)` to decrease the 
distance between  two strings, when their original distance is below some `threshold`.
The boost is equal to `min(l,  maxlength) * p * dist` where `l` denotes the 
length of their common prefix and `dist` denotes the original distance
"""
struct Winkler{S <: SemiMetric} <: SemiMetric
    dist::S
    p::Float64          # scaling factor. Default to 0.1
    threshold::Float64  # boost threshold. Default to 0.7
    maxlength::Integer      # max length of common prefix. Default to 4
    Winkler{S}(dist::S, p, threshold, maxlength) where {S <: SemiMetric} = new(dist, p, threshold, maxlength)
end

function Winkler(dist::SemiMetric; p = 0.1, threshold = 0.7, maxlength = 4)
    p * maxlength <= 1 || throw("scaling factor times maxlength of common prefix must be lower than one")
    dist = normalize(dist)
    Winkler{typeof(dist)}(dist, 0.1, 0.7, 4)
end
isnormalized(dist::Winkler) = true

function (dist::Winkler)(s1, s2, max_dist = 1.0)
    # cannot do max_dist because of boosting threshold
    score = dist.dist(s1, s2)
    if score <= 1 - dist.threshold
        l = common_prefix(s1, s2)[1]
        score -= min(l, dist.maxlength) * dist.p * score
    end
    return score
end


"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance

`TokenMax{dist}` is the minimum of the base distance `normalize(dist)`,
its [`Partial`](@ref) modifier, its [`TokenSort`](@ref) modifier, and its 
[`TokenSet`](@ref) modifier, with penalty terms depending on string lengths.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(TokenMax(RatcliffObershelp()), s1, s2)
0.05
```
"""
struct TokenMax{S <: SemiMetric} <: SemiMetric
    dist::S
    TokenMax{S}(dist::S) where {S <: SemiMetric} = new(dist)
end

function TokenMax(dist::SemiMetric)
    dist = normalize(dist)
    TokenMax{typeof(dist)}(dist)
end
isnormalized(dist::TokenMax) = true

function (dist::TokenMax)(s1::AbstractString, s2::AbstractString, max_dist = 1.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    score = dist.dist(s1, s2, max_dist)
    min_score = min(max_dist, score)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_dist = Partial(dist.dist)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
        score_partial = 1 - partial_scale * (1 - partial_dist(s1, s2, 1 - (1 - max_dist) / partial_scale))
        min_score = min(max_dist, score_partial)
        score_sort = 1 - unbase_scale * partial_scale * 
                (1 - TokenSort(partial_dist)(s1, s2, 1 - (1 - max_dist) / (unbase_scale * partial_scale)))
        max_dist = min(max_dist, score_sort)
        score_set = 1 - unbase_scale * partial_scale * 
                (1 - TokenSet(partial_dist)(s1, s2, 1 - (1 - max_dist) / (unbase_scale * partial_scale))) 
        return min(score, score_partial, score_sort, score_set)
    else
        score_sort = 1 - unbase_scale * 
                (1 - TokenSort(dist.dist)(s1, s2, 1 - (1 - max_dist) / unbase_scale))
        max_dist = min(max_dist, score_sort)
        score_set = 1 - unbase_scale * 
                (1 - TokenSet(dist.dist)(s1, s2, 1 - (1 - max_dist) / unbase_scale))
        return min(score, score_sort, score_set)
    end
end