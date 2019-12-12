
##############################################################################
##
## compare
## compare always return a value between 0 and 1.
##
##############################################################################
"""
    compare(s1::AbstractString, s2::AbstractString, dist::PreMetric)

compare returns a similarity score between the strings `s1` and `s2` based on the distance `dist`
"""
function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::Hamming; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    1.0 - evaluate(dist, s1, s2) / len2
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::Union{Jaro, RatcliffObershelp}; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    1.0 - evaluate(dist, s1, s2)
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::AbstractQGramDistance; min_score = 0.0)
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

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing},  dist::Union{Levenshtein, DamerauLevenshtein}; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    if min_score == 0.0
        return 1.0 - evaluate(dist, s1, s2) / len2
    else
        d = evaluate(dist, s1, s2; max_dist = ceil(Int, len2 * (1 - min_score)))
        out = 1.0 - d / len2
        out < min_score && return 0.0
        return out
    end
end

@deprecate compare(dist::PreMetric, s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}) compare(s1, s2, dist)

##############################################################################
##
## Winkler
##
##############################################################################
"""
   Winkler(dist::Premetric, p::Real = 0.1, boosting_threshold::Real = 0.7, l::Integer = 4)

Winkler is a `PreMetric` modifier that boosts the similarity score between two strings by a scale `p` when the strings share a common prefix with lenth lower than `l` (the boost is only applied the similarity score above `boosting_threshold`)
"""
struct Winkler{T1 <: PreMetric, T2 <: Real, T3 <: Real, T4 <: Integer} <: PreMetric
    dist::T1
    p::T2          # scaling factor. Default to 0.1
    boosting_threshold::T3      # boost threshold. Default to 0.7
    l::Integer                  # length of common prefix. Default to 4
    function Winkler(dist::T1, p::T2,  boosting_threshold::T3, l::T4) where {T1, T2, T3, T4}
        p * l >= 1 && throw("scaling factor times length of common prefix must be lower than one")
        new{T1, T2, T3, T4}(dist, p, boosting_threshold, l)
    end
end
Winkler(x) = Winkler(x, 0.1, 0.7, 4)

# hard to use min_score because of whether there is boost or not in the end
function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::Winkler; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    l = remove_prefix(s1, s2, dist.l)[1]
    # cannot do min_score because of boosting threshold
    score = compare(s1, s2, dist.dist)
    if score >= dist.boosting_threshold
        score += l * dist.p * (1 - score)
    end
    return score
end

JaroWinkler() = Winkler(Jaro(), 0.1, 0.7)

##############################################################################
##
## Partial
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
"""
   Partial(dist::Premetric)

Partial is a `PreMetric` modifier that returns the maximal similarity score between the shorter string and substrings of the longer string
"""
struct Partial{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::Partial; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 == len2 && return compare(s1, s2, dist.dist; min_score = min_score)
    len1 == 0 && return 1.0
    out = 0.0
    for x in qgram(s2, len1)
        curr = compare(s1, x, dist.dist; min_score = min_score)
        out = max(out, curr)
        min_score = max(out, min_score)
    end
    return out
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::Partial{RatcliffObershelp}; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
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

##############################################################################
##
## TokenSort
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
"""
   TokenSort(dist::Premetric)

TokenSort is a `PreMetric` modifier that adjusts for differences in word orders by reording words alphabetically.
"""
struct TokenSort{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::TokenSort; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    compare(s1, s2, dist.dist; min_score = min_score)
end

##############################################################################
##
## TokenSet
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
"""
   TokenSet(dist::Premetric)

TokenSort is a `PreMetric` modifier that adjusts for differences in word orders and word numbers by comparing the intersection of two strings with each string.
"""
struct TokenSet{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::TokenSet; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    v1 = SortedSet(split(s1))
    v2 = SortedSet(split(s2))
    v0 = intersect(v1, v2)
    s0 = join(v0, " ")
    s1 = join(v1, " ")
    s2 = join(v2, " ")
    isempty(s0) && return compare(s1, s2, dist.dist; min_score = min_score)
    dist0 = compare(s0, s1, dist.dist; min_score = min_score)
    min_score = max(min_score, dist0)
    dist1 = compare(s0, s2, dist.dist; min_score = min_score)
    min_score = max(min_score, dist1)
    dist2 = compare(s0, s2, dist.dist; min_score = min_score)
    max(dist0, dist1, dist2)
end


##############################################################################
##
## TokenMax
##
##############################################################################
"""
   TokenMax(dist::Premetric)

TokenSort is a `PreMetric` modifier that combines similarlity scores using the base distance, its Partial, TokenSort and TokenSet modifiers, with penalty terms depending on string lengths.
"""
struct TokenMax{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing}, dist::TokenMax; min_score = 0.0)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    dist0 = compare(s1, s2, dist.dist; min_score = min_score)
    min_score = max(min_score, dist0)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
        dist1 = partial_scale * compare(s1, s2, Partial(dist.dist); 
                                        min_score = min_score / partial_scale) 
        min_score = max(min_score, dist1)
        dist2 = unbase_scale * partial_scale * 
                compare(s1, s2, TokenSort(Partial(dist.dist)); 
                            min_score = min_score / (unbase_scale * partial_scale))
        min_score = max(min_score, dist2)
        dist3 = unbase_scale * partial_scale * 
                compare(s1, s2, TokenSet(Partial(dist.dist)); 
                            min_score = min_score / (unbase_scale * partial_scale)) 
        return max(dist0, dist1, dist2, dist3)
    else
        dist1 = unbase_scale * 
                compare(s1, s2, TokenSort(dist.dist); 
                            min_score = min_score / unbase_scale)
        min_score = max(min_score, dist1)
        dist2 = unbase_scale * 
                compare(s1, s2, TokenSet(dist.dist); 
                            min_score = min_score / unbase_scale) 
        return max(dist0, dist1, dist2)
    end
end