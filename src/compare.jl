\


##############################################################################
##
## compare
## compare always return a value between 0 and 1.
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    1.0 - evaluate(dist, s1, s2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::AbstractString, s2::AbstractString)
    len = max(length(s1), length(s2))
    len == 0 ? 1.0 : 1.0 - evaluate(dist, s1, s2) / len
end

function compare(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
    # When string length < q for qgram distance, returns s1 == s2
    len1 = length(s1) ; len2 = length(s2)
    min(len1, len2) <= (dist.q - 1) && return convert(Float64, s1 == s2)
    if typeof(dist) <: QGram
        1 - evaluate(dist, s1, s2) / (len1 + len2 - 2 * dist.q + 2)
    else
        1 - evaluate(dist, s1, s2)
    end
end


##############################################################################
##
## Winkler
##
##############################################################################

struct Winkler{T1 <: PreMetric, T2 <: Real, T3 <: Real} <: PreMetric
    dist::T1
    scaling_factor::T2      # scaling factor. Default to 0.1
    boosting_limit::T3      # boost threshold. Default to 0.7
end

# restrict to distance between 0 and 1
Winkler(x) = Winkler(x, 0.1, 0.7)

function compare(dist::Winkler, s1::AbstractString, s2::AbstractString)
    score = compare(dist.dist, s1, s2)
    l = common_prefix(s1, s2, 4)[1]
    # common prefix adjustment
    if score >= dist.boosting_limit
        score += l * dist.scaling_factor * (1 - score)
    end
    return score
end

##############################################################################
##
## Partial
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
struct Partial{T <: PreMetric} <: PreMetric
    dist::T
end

# general
function compare(dist::Partial, s1::AbstractString, s2::AbstractString)
    s2, len2, s1, len1 = reorder(s1, s2)
    len1 == len2 && return compare(dist.dist, s1, s2)
    len1 == 0 && return compare(dist.dist, "", "")
    iter = QGramIterator(s2, len2, len1)
    out = 0.0
    x = iterate(iter)
    while x !== nothing
        s, state = x
        curr = compare(dist.dist, s1, s)
        out = max(out, curr)
        x = iterate(iter, state)
    end
    return out
end

# Specialization for RatcliffObershelp distance
# Code follows https://github.com/seatgeek/fuzzywuzzy/blob/master/fuzzywuzzy/fuzz.py
function compare(dist::Partial{RatcliffObershelp}, s1::AbstractString, s2::AbstractString)
    s2, len2, s1, len1 = reorder(s1, s2)
    len1 == len2 && return compare(dist.dist, s1, s2)
    out = 0.0
    result = matching_blocks(s1, s2)
    for r in result
        # here I difffer from fuzz.py by making sure the substring of s2 has length len1
        s2_start = r[2] - r[1] + 1
        s2_end = s2_start + len1 - 1
        if s2_start <= 0
            s2_end += 1 - s2_start
            s2_start += 1 - s2_start
        elseif s2_end > len2
            s2_start += len2 - s2_end
            s2_end += len2 - s2_end
        end
        i2_start =  nextind(s2, 0, s2_start)
        i2_end = s2_end == len2 ? lastindex(s2) : (nextind(s2, 0, s2_end + 1) - 1)
        curr = compare(RatcliffObershelp(), s1, SubString(s2, i2_start, i2_end))
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
struct TokenSort{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(dist::TokenSort, s1::AbstractString, s2::AbstractString)
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    compare(dist.dist, s1, s2)
end

##############################################################################
##
## TokenSet
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
struct TokenSet{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(dist::TokenSet, s1::AbstractString, s2::AbstractString)
    v0, v1, v2 = _separate!(split(s1), split(s2))
    s0 = join(v0, " ")
    s1 = join(Iterators.flatten((v0, v1)), " ")
    s2 = join(Iterators.flatten((v0, v2)), " ")
    if isempty(s0)
        # otherwise compare(dist, "", "a")== 1.0 
        compare(dist.dist, s1, s2)
    else
        max(compare(dist.dist, s0, s1), 
            compare(dist.dist, s1, s2), 
            compare(dist.dist, s0, s2))        
    end
end

# separate 2 vectors in intersection, setdiff1, setdiff2 (all sorted)
function _separate!(v1::Vector, v2::Vector)
    sort!(v1)
    sort!(v2)
    out = eltype(v1)[]
    start = 1
    i1 = 0
    while i1 < length(v1)
        i1 += 1
        x = v1[i1]
        i2 = searchsortedfirst(v2, x, start, length(v2), Base.Forward)
        i2 > length(v2) && break 
        if i2 > 0 && v2[i2] == x
            deleteat!(v1, i1)
            deleteat!(v2, i2)
            push!(out, x)
            i1 -= 1
            start = i2 
        end
    end
    return out, v1, v2
end

##############################################################################
##
## TokenMax
##
##############################################################################
struct TokenMax{T <: PreMetric} <: PreMetric
    dist::T
end

function compare(dist::TokenMax, s1::AbstractString, s2::AbstractString)
    dist0 = compare(dist.dist, s1, s2)
    s2, len2, s1, len1 = reorder(s1, s2)
    unbase_scale = 0.95
    # if one string is much much shorter than the other
    if len2 >= 1.5 * len1
        # if strings are of dissimilar length, use partials
        partial = compare(Partial(dist.dist), s1, s2) 
        ptsor = compare(TokenSort(Partial(dist.dist)), s1, s2) 
        ptser = compare(TokenSet(Partial(dist.dist)), s1, s2) 
        partial_scale = len2 > (8 * len1) ? 0.6 : 0.9
        return max(dist0, 
                partial * partial_scale, 
                ptsor * unbase_scale * partial_scale, 
                ptser * unbase_scale * partial_scale)
    else
        ptsor = compare(TokenSort(dist.dist), s1, s2) 
        ptser = compare(TokenSet(dist.dist), s1, s2) 
        return max(dist0, 
                ptsor * unbase_scale, 
                ptser * unbase_scale)
    end
end