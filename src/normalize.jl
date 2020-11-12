
struct Normalized{V <: SemiMetric} <: SemiMetric
    dist::V
    max_dist::Float64
end


function (dist::Normalized{<:Hamming})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    out = dist.dist(s1, s2) / len2
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:Union{Levenshtein{Nothing}, DamerauLevenshtein{Nothing}}})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    if dist.dist isa Levenshtein
        d = Levenshtein(ceil(Int, len2 * dist.max_dist))(s1, s2)
    else
        d = DamerauLevenshtein(ceil(Int, len2 * dist.max_dist))(s1, s2)
    end
    out = d / len2
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:QGramDistance})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
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

function (dist::Normalized)(s1, s2)
    out = dist.dist(s1, s2)
    out > dist.max_dist ? 1.0 : out
end



normalize(dist::SemiMetric; max_dist = 1.0) = Normalized{typeof(dist)}(dist, max_dist)
normalize(dist::Union{Jaro, JaroWinkler}; max_dist = 1.0) = dist
normalize(dist::Partial; max_dist = 1.0) = Partial(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSort; max_dist = 1.0) = TokenSort(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSet; max_dist = 1.0) = TokenSet(normalize(dist.dist; max_dist = max_dist))
normalize(dist::Normalized; max_dist = 1.0) = Normalized{typeof(dist.dist)}(dist.dist, max_dist)

"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance

`TokenMax{dist}` normalizes the distance `dist` and returns the minimum of the distance,
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
function normalize(dist::TokenMax; max_dist = 1.0)
    dist = normalize(dist.dist; max_dist = max_dist)
    TokenMax{typeof(dist)}(dist)
end

function (dist::TokenMax)(s1::AbstractString, s2::AbstractString)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    _dist = deepcopy(dist.dist)
    max_dist = _dist.max_dist
    score = _dist(s1, s2)
    min_score = min(max_dist, score)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
       _dist = Normalized(_dist.dist, 1 - (1 - max_dist) / partial_scale)
        score_partial = 1 - partial_scale * (1 - Partial(_dist)(s1, s2))
        min_score = min(max_dist, score_partial)
       _dist = Normalized(_dist.dist, 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_sort = 1 - unbase_scale * partial_scale * (1 - TokenSort(Partial(_dist))(s1, s2))
        max_dist = min(max_dist, score_sort)
       _dist = Normalized(_dist.dist, 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_set = 1 - unbase_scale * partial_scale * (1 - TokenSet(Partial(_dist))(s1, s2)) 
        out = min(score, score_partial, score_sort, score_set)
    else
       _dist = Normalized(_dist.dist, 1 - (1 - max_dist) / unbase_scale)
        score_sort = 1 - unbase_scale * (1 - TokenSort(_dist)(s1, s2))
        max_dist = min(max_dist, score_sort)
       _dist = Normalized(_dist.dist,  1 - (1 - max_dist) / unbase_scale)
        score_set = 1 - unbase_scale * (1 - TokenSet(_dist)(s1, s2))
        out = min(score, score_sort, score_set)
    end
    out > max_dist ? 1.0 : out
end
