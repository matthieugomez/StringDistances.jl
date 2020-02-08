"""
    compare(s1::AbstractString, s2::AbstractString, dist::StringDistance)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the string distance `dist`.

### Examples
```julia-repl
julia> compare("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function compare(s1, s2, dist::Union{Jaro, RatcliffObershelp}; min_score = 0.0)
    1.0 - evaluate(dist, s1, s2)
end

function compare(s1, s2,  dist::Union{Levenshtein, DamerauLevenshtein}; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    d = evaluate(dist, s1, s2; max_dist = ceil(Int, len2 * (1 - min_score)))
    out = 1.0 - d / len2
    out < min_score ? 0.0 : out
end

function compare(s1, s2, dist::QGramDistance; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    # When string length < q for qgram distance, returns s1 == s2
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 <= dist.q - 1 && return convert(Float64, s1 == s2)
    if typeof(dist) <: QGram
        1.0 - evaluate(dist, s1, s2) / (len1 + len2 - 2 * dist.q + 2)
    else
        1.0 - evaluate(dist, s1, s2)
    end
end

"""
   Winkler(dist::StringDistance; p::Real = 0.1, threshold::Real = 0.7, maxlength::Integer = 4)

Creates the `Winkler{dist, p, threshold, maxlength}` distance

`Winkler{dist, p, threshold, length)` modifies the string distance `dist` to boost the 
similarity score between  two strings, when their original similarity score is above some `threshold`.
The boost is equal to `min(l,  maxlength) * p * (1 - score)` where `l` denotes the 
length of their common prefix and `score` denotes the original score
"""
struct Winkler{S <: SemiMetric} <: SemiMetric
    dist::S
    p::Float64          # scaling factor. Default to 0.1
    threshold::Float64  # boost threshold. Default to 0.7
    maxlength::Integer      # max length of common prefix. Default to 4
end

function Winkler(dist::StringDistance; p = 0.1, threshold = 0.7, maxlength = 4)
    p * maxlength <= 1 || throw("scaling factor times maxlength of common prefix must be lower than one")
    Winkler(dist, 0.1, 0.7, 4)
end

function compare(s1, s2, dist::Winkler; min_score = 0.0)
    # cannot do min_score because of boosting threshold
    score = compare(s1, s2, dist.dist)
    if score >= dist.threshold
        l = common_prefix(s1, s2)[1]
        score += min(l, dist.maxlength) * dist.p * (1 - score)
    end
    return score
end


"""
   Partial(dist::StringDistance)

Creates the `Partial{dist}` distance

`Partial{dist}` modifies the string distance `dist` to return the 
maximal similarity score  between the shorter string and substrings of the longer string

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> compare(s1, s2, Partial(RatcliffObershelp()))
0.4516129032258065
```
"""
struct Partial{S <: SemiMetric} <: SemiMetric
    dist::S
end

function compare(s1, s2, dist::Partial; min_score = 0.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return compare(s1, s2, dist.dist; min_score = min_score)
    len1 == 0 && return 1.0
    out = 0.0
    for x in qgrams(s2, len1)
        curr = compare(s1, x, dist.dist; min_score = min_score)
        out = max(out, curr)
        min_score = max(out, min_score)
    end
    return out
end

function compare(s1::AbstractString, s2::AbstractString, dist::Partial{RatcliffObershelp}; min_score = 0.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return compare(s1, s2, dist.dist)
    out = 0.0
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
        i2_start = nextind(s2, 0, s2_start)
        i2_end = nextind(s2, 0, s2_end)
        curr = compare(s1, SubString(s2, i2_start, i2_end), RatcliffObershelp())
        out = max(out, curr)
    end
    return out
end

"""
   TokenSort(dist::StringDistance)

Creates the `TokenSort{dist}` distance

`TokenSort{dist}` modifies the string distance `dist` to adjust for differences 
in word orders by reording words alphabetically.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> compare(s1, s2, TokenSort(RatcliffObershelp()))
1.0
```
"""
struct TokenSort{S <: SemiMetric} <: SemiMetric
    dist::S
end

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function compare(s1, s2, dist::TokenSort; min_score = 0.0)
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    compare(s1, s2, dist.dist; min_score = min_score)
end


"""
   TokenSet(dist::StringDistance)

Creates the `TokenSet{dist}` distance

`TokenSet{dist}` modifies the string distance `dist` to adjust for differences 
in word orders and word numbers by comparing the intersection of two strings with each string.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> compare(s1, s2, TokenSet(RatcliffObershelp()))
1.0
```
"""
struct TokenSet{S <: SemiMetric} <: SemiMetric
    dist::S
end

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function compare(s1, s2, dist::TokenSet; min_score = 0.0)
    v1 = unique!(sort!(split(s1)))
    v2 = unique!(sort!(split(s2)))
    v0 = intersect(v1, v2)
    s0 = join(v0, " ")
    s1 = join(v1, " ")
    s2 = join(v2, " ")
    isempty(s0) && return compare(s1, s2, dist.dist; min_score = min_score)
    score_01 = compare(s0, s1, dist.dist; min_score = min_score)
    min_score = max(min_score, score_01)
    score_02 = compare(s0, s2, dist.dist; min_score = min_score)
    min_score = max(min_score, score_02)
    score_12 = compare(s1, s2, dist.dist; min_score = min_score)
    max(score_01, score_02, score_12)
end


"""
   TokenMax(dist::StringDistance)

Creates the `TokenMax{dist}` distance

`TokenMax{dist}` combines similarity scores of the base distance `dist`,
its [`Partial`](@ref) modifier, its [`TokenSort`](@ref) modifier, and its 
[`TokenSet`](@ref) modifier, with penalty terms depending on string lengths.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> compare(s1, s2, TokenMax(RatcliffObershelp()))
0.95
```
"""
struct TokenMax{S <: SemiMetric} <: SemiMetric
    dist::S
end

function compare(s1, s2, dist::TokenMax; min_score = 0.0)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    score = compare(s1, s2, dist.dist; min_score = min_score)
    min_score = max(min_score, score)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
        score_partial = partial_scale * compare(s1, s2, Partial(dist.dist); 
                                        min_score = min_score / partial_scale) 
        min_score = max(min_score, score_partial)
        score_sort = unbase_scale * partial_scale * 
                compare(s1, s2, TokenSort(Partial(dist.dist)); 
                            min_score = min_score / (unbase_scale * partial_scale))
        min_score = max(min_score, score_sort)
        score_set = unbase_scale * partial_scale * 
                compare(s1, s2, TokenSet(Partial(dist.dist)); 
                            min_score = min_score / (unbase_scale * partial_scale)) 
        return max(score, score_partial, score_sort, score_set)
    else
        score_sort = unbase_scale * 
                compare(s1, s2, TokenSort(dist.dist); 
                            min_score = min_score / unbase_scale)
        min_score = max(min_score, score_sort)
        score_set = unbase_scale * 
                compare(s1, s2, TokenSet(dist.dist); 
                            min_score = min_score / unbase_scale) 
        return max(score, score_sort, score_set)
    end
end