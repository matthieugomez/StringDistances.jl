"""
   Partial(dist)

Creates the `Partial{dist}` distance.

`Partial{dist}`  returns the  minimum distance  between the shorter string and substrings of the longer string that have a length equal to the shorter string.

See http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> Partial(RatcliffObershelp())(s1, s2)
0.5483870967741935
```
"""
struct Partial{S <: Union{StringSemiMetric, StringMetric}} <: StringSemiMetric
    dist::S
end

function (dist::Partial)(s1, s2; max_dist = nothing)
    (s1 === missing) | (s2 === missing) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    out = dist.dist(s1, s2; max_dist = max_dist)
    max_dist0 = (max_dist !== nothing) ? min(max_dist, out) : out
    ((len1 == 0) | (len1 == len2)) && return out
    for x in qgrams(s2, len1)
        curr = dist.dist(s1, x; max_dist = max_dist0)
        out = min(out, curr)
        max_dist0 = min(max_dist0, curr)
    end
    return out
end

# specialized (faster) version for RatcliffObershelp
function (dist::Partial{<: Union{RatcliffObershelp, Normalized{RatcliffObershelp}}})(s1, s2; max_dist = nothing)
    (s1 === missing) | (s2 === missing) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return dist.dist(s1, s2)
    out = 1.0
    for s2_start in matching_blocks(s1, s2, 1, 1, len1, len2)
        # Make sure the substring of s2 has length len1
        if s2_start < 1
            s2_start = 1
        elseif s2_start + len1 - 1 > len2
            s2_start += len2 - (s2_start + len1 - 1)
        end
        n_matched = length_matching_blocks(s1, s2, 1, s2_start, len1, s2_start + len1 - 1)
        curr = 1 - 2 * n_matched / (len1 + len1)
        out = min(out, curr)
    end
    return out
end

function matching_blocks(s1, s2, start1::Integer, start2::Integer, end1::Integer, end2::Integer)
    x = Set{Int}()
    p = zeros(Int, max(end1 - start1, end2 - start2) + 1)
    matching_blocks!(x, p, s1, s2, start1, start2, end1, end2)
end

function matching_blocks!(x::Set{Int}, p::Vector{Int}, s1, s2, start1::Integer, start2::Integer, end1::Integer, end2::Integer)
    j1, j2, len = longest_common_pattern!(p, s1, s2, start1, start2, end1, end2)
    len == 0 && return x
    push!(x, j2 - j1 + 1)
    matching_blocks!(x, p, s1, s2, start1, start2, j1 - 1, j2 - 1)
    matching_blocks!(x, p, s1, s2, j1 + len, j2 + len, end1, end2)
    return x
end

Normalized(dist::Partial) = Normalized{typeof(Partial(Normalized(dist.dist)))}(Partial(Normalized(dist.dist)))
"""
   TokenSort(dist)

Creates the `TokenSort{dist}` distance.

`TokenSort{dist}` returns the distance between strings after reording words alphabetically.
See http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

It is only defined on AbstractStrings.

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s1 = "New York Mets vs Atlanta Braves"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> TokenSort(RatcliffObershelp())(s1, s2)
0.0
```
"""
struct TokenSort{S <: Union{StringSemiMetric, StringMetric}} <: StringSemiMetric
    dist::S
end

function (dist::TokenSort)(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}; max_dist = nothing)
    (s1 === missing) | (s2 === missing) && return missing
    f = s -> join(sort!(split(s)), " ")
    dist.dist(f(s1), f(s2); max_dist = max_dist)
end

Normalized(dist::TokenSort) = Normalized{typeof(TokenSort(Normalized(dist.dist)))}(TokenSort(Normalized(dist.dist)))
"""
   TokenSet(dist)

Creates the `TokenSet{dist}` distance, which is only defined on AbstractStrings.

`TokenSet{dist}` returns the minimum the distances between:
[SORTED_INTERSECTION]
[SORTED_INTERSECTION] + [SORTED_REST_OF_STRING1]
[SORTED_INTERSECTION] + [SORTED_REST_OF_STRING2]
See: http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> TokenSet(RatcliffObershelp())(s1, s2)
0.0
```
"""
struct TokenSet{S <: Union{StringSemiMetric, StringMetric}} <: StringSemiMetric
    dist::S
end

function (dist::TokenSet)(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}; max_dist = nothing)
    (s1 === missing) | (s2 === missing) && return missing
    v1 = unique!(sort!(split(s1)))
    v2 = unique!(sort!(split(s2)))
    v0 = intersect(v1, v2)
    s0 = join(v0, " ")
    s1 = join(v1, " ")
    s2 = join(v2, " ")
    isempty(s0) && return dist.dist(s1, s2; max_dist = max_dist)
    min(dist.dist(s0, s1; max_dist = max_dist),
        dist.dist(s0, s2; max_dist = max_dist),
        dist.dist(s1, s2; max_dist = max_dist))
end

Normalized(dist::TokenSet) = Normalized{typeof(TokenSet(Normalized(dist.dist)))}(TokenSet(Normalized(dist.dist)))
"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance, which is only defined on AbstractStrings.

`TokenMax{dist}` normalizes the distance `dist` and returns the minimum of the distance,
its [`Partial`](@ref) modifier, its [`TokenSort`](@ref) modifier, and its 
[`TokenSet`](@ref) modifier, with penalty terms depending on the strings lengths.


### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(TokenMax(RatcliffObershelp()), s1, s2)
0.05
```
"""
struct TokenMax{S <: Normalized} <: StringSemiMetric
    dist::S
end
TokenMax(dist::Union{StringSemiMetric, StringMetric}) = TokenMax(Normalized(dist))

function (dist::TokenMax)(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}; max_dist = 1.0)
    (s1 === missing) | (s2 === missing) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    dist0 = dist.dist
    out = dist0(s1, s2; max_dist = max_dist)
    max_dist = min(max_dist, out)
    scale = 0.95
    # if one string is much shorter than the other, use partial
    if len2 >= 1.5 * len1
        dist0 = Partial(dist0)
        pscale = 0.9
        pout = 1 - pscale *  (1 - dist0(s1, s2; max_dist = 1 - (1 - max_dist) / pscale))
        out = min(out, pout)
        max_dist = min(max_dist, pout)
        scale *= pscale
    end
    out_sort = 1 - scale * (1 - TokenSort(dist0)(s1, s2; max_dist = 1 - (1 - max_dist) / scale))
    max_dist = min(max_dist, out_sort)
    out_set = 1 - scale * (1 - TokenSet(dist0)(s1, s2; max_dist = 1 - (1 - max_dist) / scale))
    out = min(out, out_sort, out_set)
    out > max_dist ? 1.0 : out
end
