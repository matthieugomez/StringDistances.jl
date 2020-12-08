"""
    update_max_dist(dist, max_dist) -> StringDistance

Given a `StringDistance` `dist` which supports the `max_dist` field, return
a similar object with the field updated to the value passed in. If the distance does
not support this field, `dist` is returned.
"""
update_max_dist

# Why not use `update_max_dist(::T, max_dist) where T <: Union{DamerauLevenshtein, Levenshtein, Hamming}`
# instead of `@eval`? Because then `T` will be e.g. `DamerauLevenshtein{Nothing}`, and we cannot then
# do `T(1)` to create a new one.
for T in (:DamerauLevenshtein, :Levenshtein, :Hamming)
    @eval function update_max_dist(::$T, max_dist)
        return $T(max_dist)
    end
end

for T in (:Normalized, :TokenMax)
    @eval function update_max_dist(dist::$T, max_dist)
        return $T(dist.dist, max_dist)
    end
end

for T in (:Partial, :TokenSort, :TokenSet)
    @eval function update_max_dist(dist::$T, max_dist)
        return $T(update_max_dist(dist.dist, max_dist))
    end
end
update_max_dist(d::Any, max_dist) = d

"""
    get_max_dist(dist) -> Number

Given a `StringDistance` `dist` which supports the `max_dist` field, return
the value of the field. If the object does not support `max_dist`, then
return nothing.
"""
get_max_dist

get_max_dist(dist::Union{Levenshtein, DamerauLevenshtein, Hamming, Normalized, TokenMax}) = dist.max_dist
get_max_dist(dist::Union{Partial, TokenSort, TokenSet}) = get_max_dist(dist.dist)
get_max_dist(::Any) = nothing

"""
    findnearest_partial(needle, haystack, dist) -> (d, inds)

`Partial(dist)(needle, haystack)` returns
the closest distance `d` between `needle` and any segment of `haystack` (of equal length to that of `needle`). `findnearest_partial` returns the same value, but also
returns the first set of indices at which an optimal partial match was found. If `dist` supports a `max_dist`
field, and no match was found with distance at most `max_dist`, then returns an empty range (`1:0`) for
the indices. Requires `haystack` to be indexable (e.g. an `AbstractString`).

See also [`Partial`](@ref) and [`findall_partial`](@ref).
"""
findnearest_partial

# unwrap `Partial`s since we compare as partials anyway.
findnearest_partial(s1, s2, dist::Partial) = findnearest_partial(s1, s2, dist.dist)

function findnearest_partial(s1, s2, dist)
    max_dist = get_max_dist(dist)
    s1, s2 = enforce_shorter_first(s1, s2)  

    if max_dist === nothing
        # return something larger than any possible distance,
        # but not e.g. `typemax(Int)` which will lead to overflows,
        # and with an integer type, since we need to be able to
        # construct e.g. `Levenshtein` distances with this parameter.
        max_dist = 10*length(s2)
    end

    len1, len2 = length(s1), length(s2)
    len1 == len2 && return dist(s1, s2), firstindex(s2):lastindex(s2)
    out = max_dist + 1
    len1 == 0 && return out, 1:0
    out_idx = 0
    for (i, x) in enumerate(qgrams(s2, len1))
        curr = dist(s1, x)
        out_idx = ifelse(curr < out, i, out_idx)
        out = min(out, curr)
        max_dist = max_dist === nothing ? out : min(out, max_dist)
        dist = update_max_dist(dist, max_dist)
    end

    if out_idx == 0
        # return more obvious invalid range if a match isn't found without exceeding `max_dist`
        return out, 1:0
    else
        return out, _slice_inds(s2, out_idx, out_idx + len1 - 1)
    end
end


"""
    findall_partial(needle, haystack, dist; max_dist = StringDistances.get_max_dist(dist)) -> Vector{Tuple{T,UnitRange}}

Searches for occurrences of `needle` in `haystack` that differ by at most `max_dist` according to the distance measure `dist`. Only considers sequential segments of `haystack` of equal length to that of `needle`.

Returns a vector of tuples, each corresponding to a match found. The first entry gives the distance of the match, and the second entry gives the indices of `haystack` corresponding to the match. Matches may overlap. Requires `haystack` to be indexable (e.g. an `AbstractString`).

See also [`Partial`](@ref) and [`findnearest_partial`](@ref).
"""
findall_partial

findall_partial(s1, s2, dist::Partial; max_dist = get_max_dist(dist)) = findall_partial(s1, s2, dist.dist; max_dist = max_dist)

function findall_partial(s1, s2, dist; max_dist = get_max_dist(dist))
    if max_dist === nothing
        throw(ArgumentError("`dist` does not have a `max_dist` set and one was not passed to `findall_partial`."))
    end

    s1, s2 = enforce_shorter_first(s1, s2)  
    T = Distances.result_type(dist, s1, s2)

    dist = update_max_dist(dist, max_dist)
    len1, len2 = length(s1), length(s2)
    matches = Tuple{T,UnitRange}[]
    len1 == 0 && return matches

    if len1 == len2
        curr = dist(s1, s2)
        if curr <= max_dist
            push!(matches, (curr, firstindex(s2):lastindex(s2)))
        end
        return matches
    end

    for (i, x) in enumerate(qgrams(s2, len1))
        curr = dist(s1, x)
        if curr <= max_dist
            inds = _slice_inds(s2, i, i + len1 - 1)
            push!(matches, (curr, inds))
        end
    end
    return matches
end
