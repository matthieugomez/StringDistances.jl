struct Normalized{V <: Union{StringSemiMetric, StringMetric}} <: StringSemiMetric
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
normalize(dist::Union{StringSemiMetric, StringMetric}; max_dist = 1.0) = Normalized{typeof(dist)}(dist, max_dist)
normalize(dist::Normalized; max_dist = 1.0) = Normalized(dist.dist, max_dist)


"""
    compare(s1, s2, dist)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the distance `dist`.

### Examples
```julia-repl
julia> compare("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function compare(s1, s2, dist::Union{StringSemiMetric, StringMetric}; min_score = 0.0)
    1 - normalize(dist, max_dist = 1 - min_score)(s1, s2)
end 
