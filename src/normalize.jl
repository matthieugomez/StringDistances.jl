# Normalized is basically wrapper like Symmetric 

"""
   Normalized(dist::Union{StringSemiMetric, StringMetric})

Creates a normalized distance. The normalized distance always return a Float64 between 0.0 and 1.0 (or a missing if one of the argument is missing)

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> Levenshtein()(s1, s2)
25
julia> StringDistances.Normalized(Levenshtein())(s1, s2)
0.8064 
```
"""
struct Normalized{T <: Union{StringSemiMetric, StringMetric}} <: NormalizedStringSemiMetric
    dist::T
end
Normalized(dist::Union{StringSemiMetric, StringMetric}) = Normalized{typeof(dist)}(dist)
Normalized(dist::Normalized) = dist


# this basically says that all distances are considered to be normalized by default
function (dist::Normalized)(s1, s2; max_dist = 1.0)
    dist.dist(s1, s2; max_dist = max_dist)
end



function (dist::Normalized{<:Union{Hamming, DamerauLevenshtein}})(s1, s2; max_dist = 1.0)
    (s1 === missing) | (s2 === missing) && return missing
    isempty(s1) && isempty(s2) && return 0.0
    out = dist.dist(s1, s2) / length(s2)
    max_dist !== nothing && out > max_dist && return 1.0
    return out
end

function (dist::Normalized{<:Union{Levenshtein, OptimalStringAlignement}})(s1, s2; max_dist = 1.0)
    (s1 === missing) | (s2 === missing) && return missing
    isempty(s1) && isempty(s2) && return 0.0
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    if max_dist === nothing || max_dist == 1.0
        d = dist.dist(s1, s2)
    else
        d = dist.dist(s1, s2; max_dist = ceil(Int, len2 * max_dist))
    end
    out = d / len2
    max_dist !== nothing && out > max_dist && return 1.0
    return out
end

function (dist::Normalized{<:AbstractQGramDistance})(s1, s2; max_dist = 1.0)
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
    max_dist !== nothing && out > max_dist && return 1.0
    return out
end