struct Normalize{S <: SemiMetric} <: SemiMetric
    dist::S
end
"""
   normalize(dist::SemiMetric)

   Normalize a metric, so that `evaluate` always return a Float64 between 0 and 1 (or a `missing` if one element is missing)
"""
# also a normalized distance always accept a third argument, max_dist.

normalize(dist::SemiMetric) = Normalize{typeof(dist)}(dist)

function (dist::Normalize{<: Union{Levenshtein, DamerauLevenshtein}})(s1, s2, max_dist = 1.0)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    d = dist.dist(s1, s2, ceil(Int, len2 * max_dist))
    out = d / len2
    out > max_dist ? 1.0 : out
end

function (dist::Normalize{<: QGramDistance})(s1, s2, max_dist = 1.0)
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

function (dist::Normalize)(s1, s2, max_dist = 1.0)
    dist.dist(s1, s2)
end

"""
   Winkler(dist; p::Real = 0.1, threshold::Real = 0.7, maxlength::Integer = 4)

Creates the `Winkler{dist, p, threshold, maxlength}` distance.

`Winkler{dist, p, threshold, length)` modifies the string distance `dist` to decrease the 
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
    Winkler{typeof(normalize(dist))}(normalize(dist), 0.1, 0.7, 4)
end
normalize(dist::Winkler) = dist

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
   Partial(dist)

Creates the `Partial{dist}` distance.

`Partial{dist}` modifies the string distance `dist` to return the 
minimum distance  between the shorter string and substrings of the longer string

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(Partial(RatcliffObershelp()), s1, s2)
0.5483870967741935
```
"""
struct Partial{S <: SemiMetric} <: SemiMetric
    dist::S
    Partial{S}(dist::S) where {S <: SemiMetric} = new(dist)
end
Partial(dist::SemiMetric) = Partial{typeof(normalize(dist))}(normalize(dist))
normalize(dist::Partial) = dist

function (dist::Partial)(s1, s2, max_dist = 1.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return dist.dist(s1, s2, max_dist)
    len1 == 0 && return 1.0
    out = 1.0
    for x in qgrams(s2, len1)
        curr = dist.dist(s1, x, max_dist)
        out = min(out, curr)
        max_dist = min(out, max_dist)
    end
    return out
end

function (dist::Partial{Normalize{RatcliffObershelp}})(s1, s2, max_dist = 1.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return dist.dist(s1, s2)
    out = 1.0
    for r in matching_blocks(s1, s2)
        # Make sure the substring of s2 has length len1
        s2_start = r[2] - r[1] + 1
        s2_end = s2_start + len1 - 1
        if s2_start <= 0
            s2_end += 1 - s2_start
            s2_start += 1 - s2_start
        elseif s2_end > len2
            s2_start += len2 - s2_end
            s2_end += len2 - s2_end
        end
        curr = dist.dist(s1, _slice(s2, s2_start - 1, s2_end))
        out = min(out, curr)
    end
    return out
end

"""
   TokenSort(dist)

Creates the `TokenSort{dist}` distance.

`TokenSort{dist}` modifies the string distance `dist` to adjust for differences 
in word orders by reording words alphabetically.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(TokenSort(RatcliffObershelp()), s1, s2)
0.0
```
"""
struct TokenSort{S <: SemiMetric} <: SemiMetric
    dist::S
    TokenSort{S}(dist::S) where {S <: SemiMetric} = new(dist)
end
TokenSort(dist::SemiMetric) = TokenSort{typeof(normalize(dist))}(normalize(dist))
normalize(dist::TokenSort) = dist

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function (dist::TokenSort)(s1::AbstractString, s2::AbstractString, max_dist = 1.0)
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    dist.dist(s1, s2, max_dist)
end


"""
   TokenSet(dist)

Creates the `TokenSet{dist}` distance.

`TokenSet{dist}` modifies the string distance `dist` to adjust for differences 
in word orders and word numbers by comparing the intersection of two strings with each string.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(TokenSet(RatcliffObershelp()), s1, s2)
0.0
```
"""
struct TokenSet{S <: SemiMetric} <: SemiMetric
    dist::S
    TokenSet{S}(dist::S) where {S <: SemiMetric} = new(dist)
end
TokenSet(dist::SemiMetric) = TokenSet{typeof(normalize(dist))}(normalize(dist))
normalize(dist::TokenSet) = dist

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function (dist::TokenSet)(s1::AbstractString, s2::AbstractString, max_dist = 1.0)
    v1 = unique!(sort!(split(s1)))
    v2 = unique!(sort!(split(s2)))
    v0 = intersect(v1, v2)
    s0 = join(v0, " ")
    s1 = join(v1, " ")
    s2 = join(v2, " ")
    isempty(s0) && return dist.dist(s1, s2, max_dist)
    score_01 = dist.dist(s0, s1, max_dist)
    max_dist = min(max_dist, score_01)
    score_02 = dist.dist(s0, s2, max_dist)
    max_dist = min(max_dist, score_02)
    score_12 = dist.dist(s1, s2, max_dist)
    min(score_01, score_02, score_12)
end


"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance

`TokenMax{dist}` is the minimum of the base distance `dist`,
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

TokenMax(dist::SemiMetric) = TokenMax{typeof(normalize(dist))}(normalize(dist))
normalize(dist::TokenMax) = dist

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