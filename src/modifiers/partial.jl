##############################################################################
##
## Partial
## From the Python module fuzzywuzzy
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
type Partial{T <: PreMetric} <: PreMetric
    dist::T
end

# general
function compare(dist::Partial, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    len1 == len2 && return compare(dist.dist, s1, s2, len1, len2)
    len1 == 0 && return compare(dist.dist, "", "", 0, 0)
    iter = QGramIterator(s2, len2, len1)
    state = start(iter)
    s, state = next(iter, state)
    out = compare(dist.dist, s1, s)
    while !done(iter, state)
        s, state = next(iter, state)
        curr = compare(dist.dist, s1, s)
        out = max(out, curr)
    end
    return out
end

# Specialization for RatcliffObershelp distance
# Code: https://github.com/seatgeek/fuzzywuzzy/blob/master/fuzzywuzzy/fuzz.py
function compare(dist::Partial{RatcliffObershelp}, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    len1 == len2 && return compare(dist.dist, s1, s2, len1, len2)
    out = 0.0
    result = matching_blocks(s1, s2)
    for r in result
        s2_start = max(1, r[2] - r[1] + 1)
        s2_end = s2_start + len1 - 1
        i2_start =  chr2ind(s2, s2_start)
        i2_end = s2_end == len2 ? endof(s2) : (chr2ind(s2, s2_end + 1) - 1)
        curr = compare(RatcliffObershelp(), s1, SubString(s2, i2_start, i2_end), len1, len1)
        out = max(out, curr)
    end
    return out
end
