"""
   Partial(dist)

Creates the `Partial{dist}` distance.

`Partial{dist}`  returns the  minimum distance  between the shorter string and substrings of the longer string (of the size of the shorter stirng)

See: http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

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
end

function (dist::Partial)(s1, s2)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    out = dist.dist(s1, s2)
    ((len1 == 0) | (len1 == len2)) && return out
    for x in qgrams(s2, len1)
        curr = dist.dist(s1, x)
        out = min(out, curr)
    end
    return out
end

function (dist::Partial{RatcliffObershelp})(s1, s2)
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return dist.dist(s1, s2)
    out = 1.0
    for r in matching_blocks(s1, s2)
        # Make sure the substring of s2 has length len1
        s2_start = r[2] - r[1] + 1
        s2_end = s2_start + len1 - 1
        if s2_start < 1
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

`TokenSort{dist}` returns the distance between strings after reording words alphabetically.

See: http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

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
end

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function (dist::TokenSort)(s1::AbstractString, s2::AbstractString)
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    out = dist.dist(s1, s2)
end

"""
   TokenSet(dist)

Creates the `TokenSet{dist}` distance.

`TokenSet{dist}` compares the intersection of two strings with each string, after reording words alphabetically

See: http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

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
end

# http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
function (dist::TokenSet)(s1::AbstractString, s2::AbstractString)
    v1 = unique!(sort!(split(s1)))
    v2 = unique!(sort!(split(s2)))
    v0 = intersect(v1, v2)
    s0 = join(v0, " ")
    s1 = join(v1, " ")
    s2 = join(v2, " ")
    isempty(s0) && return dist.dist(s1, s2)
    score_01 = dist.dist(s0, s1)
    score_02 = dist.dist(s0, s2)
    score_12 = dist.dist(s1, s2)
    min(score_01, score_02, score_12)
end


