struct Normalized{V <: SemiMetric} <: SemiMetric
    dist::V
    max_dist::Float64
end

function (dist::Normalized{<: Union{Jaro, JaroWinkler, RatcliffObershelp}})(s1, s2)
    out = dist.dist(s1, s2)
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:Union{Hamming, DamerauLevenshtein}})(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    isempty(s1) && isempty(s2) && return 0.0
    out = dist.dist(s1, s2) / length(s2)
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:Union{Levenshtein{Nothing}, OptimalStringAlignement{Nothing}}})(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    isempty(s1) && isempty(s2) && return 0.0
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    if dist.dist isa Levenshtein
        d = Levenshtein(ceil(Int, len2 * dist.max_dist))(s1, s2)
    else
        d = OptimalStringAlignement(ceil(Int, len2 * dist.max_dist))(s1, s2)
    end
    out = d / len2
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:AbstractQGramDistance})(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    # When string length < q for qgram distance, returns s1 == s2
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 <= dist.dist.q - 1 && return convert(Float64, s1 != s2)
    if dist.dist isa QGram
        out = dist.dist(s1, s2) / (len1 + len2 - 2 * dist.dist.q + 2)
    else
        out = dist.dist(s1, s2)
    end
    out > dist.max_dist ? 1.0 : out
end


"""
   normalize(dist)

Creates a normalized distance. The distance always return a Float64 between 0.0 and 1.0 (or a missing if one of the argument is missing)

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> Levenshtein()(s1, s2)
25
julia> StringDistances.normalize(Levenshtein())(s1, s2)
0.8064 
```
"""
normalize(dist::SemiMetric; max_dist = 1.0) = Normalized{typeof(dist)}(dist, max_dist)
normalize(dist::Partial; max_dist = 1.0) = Partial(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSort; max_dist = 1.0) = TokenSort(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSet; max_dist = 1.0) = TokenSet(normalize(dist.dist; max_dist = max_dist))
normalize(dist::Normalized; max_dist = 1.0) = Normalized(dist.dist, max_dist)


"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance.

`TokenMax{dist}` normalizes the distance `dist` and returns the minimum of the distance,
its [`Partial`](@ref) modifier, its [`TokenSort`](@ref) modifier, and its 
[`TokenSet`](@ref) modifier, with penalty terms depending on the iterator length.


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
    max_dist::Float64
end
TokenMax(dist::SemiMetric; max_dist = 1.0) = TokenMax(dist, max_dist)
normalize(dist::TokenMax; max_dist = 1.0) = TokenMax(dist.dist, max_dist)

function (dist::TokenMax)(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    max_dist = dist.max_dist
    dist0 = normalize(dist.dist; max_dist = max_dist)
    score = dist0(s1, s2)
    min_score = min(max_dist, score)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / partial_scale)
        score_partial = 1 - partial_scale * (1 - Partial(dist0)(s1, s2))
        min_score = min(max_dist, score_partial)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_sort = 1 - unbase_scale * partial_scale * (1 - TokenSort(Partial(dist0))(s1, s2))
        max_dist = min(max_dist, score_sort)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_set = 1 - unbase_scale * partial_scale * (1 - TokenSet(Partial(dist0))(s1, s2)) 
        out = min(score, score_partial, score_sort, score_set)
    else
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / unbase_scale)
        score_sort = 1 - unbase_scale * (1 - TokenSort(dist0)(s1, s2))
        max_dist = min(max_dist, score_sort)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / unbase_scale)
        score_set = 1 - unbase_scale * (1 - TokenSet(dist0)(s1, s2))
        out = min(score, score_sort, score_set)
    end
    out > max_dist ? 1.0 : out
end
