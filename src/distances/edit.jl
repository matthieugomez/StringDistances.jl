"""
    Hamming()

Creates the Hamming distance

The Hamming distance is defined as the number of characters that do not match
"""
struct Hamming{V <: Union{Int, Nothing}} <: SemiMetric
   max_dist::V
end
Hamming() = Hamming(nothing)

function (dist::Hamming{T})(s1, s2) where {T}
    (s1 === missing) | (s2 === missing) && return missing
    out = abs(length(s2) - length(s1))
    for (ch1, ch2) in zip(s1, s2)
        out += ch1 != ch2
        if T <: Int
            out > dist.max_dist && return dist.max_dist + 1
        end
    end
    return out
end


"""
    Jaro()

Creates the Jaro distance

The Jaro distance is defined as


``1 - (m / |s1| + m / |s2| + (m - t) / m) / 3``

where ``m`` is the number of matching characters and 
``t`` is half the number of transpositions.
"""
struct Jaro <: SemiMetric end

## http://alias-i.com/lingpipe/docs/api/com/aliasi/spell/JaroWinklerDistance.html
function (dist::Jaro)(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    end
    # If both iterators empty, formula in Wikipedia gives 1, but it makes more sense to set it to s1 == s2
    len2 > 0 || return Float64(s1 == s2)
    d = max(0, div(len2, 2) - 1)
    flag = fill(false, len2)
    ch1_match = Vector{eltype(s1)}()
    for (i1, ch1) in enumerate(s1)
        for (i2, ch2) in enumerate(s2)
            # for each character in s1, greedy search of matching character in s2 within a distance d
            i2 >= i1 - d || continue
            i2 <= i1 + d || break
            if ch1 == ch2 && !flag[i2] 
                flag[i2] = true
                push!(ch1_match, ch1)
                break
            end
        end
    end
    if isempty(ch1_match)
        return 1.0
    else
        #  m counts number matching characters
        m = length(ch1_match)
        # t/2 counts number transpositions
        t = 0
        i1 = 0
        for (i2, ch2) in enumerate(s2)
            if flag[i2]
                i1 += 1
                @inbounds t += ch2 != ch1_match[i1]
            end
        end
        return 1 - (m / len1 + m / len2 + (m - 0.5 * t) / m) / 3
    end
end


"""
    JaroWinkler(;p = 0.1, threshold = 0.3, maxlength = 4)

Creates the JaroWinkler distance

The JaroWinkler distance is defined as the Jaro distance, which is multiplied by
``(1-min(l,  maxlength) * p)`` as long as it is lower than `threshold`, and where `l` denotes the length of the common prefix.
"""
struct JaroWinkler <: SemiMetric
    p::Float64          # scaling factor. Default to 0.1
    threshold::Float64  # boost limit. Default to 0.3
    maxlength::Integer  # max length of common prefix. Default to 4
end

JaroWinkler(; p = 0.1, threshold = 0.3, maxlength = 4) = JaroWinkler(p, threshold, maxlength)

## http://alias-i.com/lingpipe/docs/api/com/aliasi/spell/JaroWinklerDistance.html
function (dist::JaroWinkler)(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    out = Jaro()(s1, s2)
    if out <= dist.threshold
        l = common_prefix(s1, s2)[1]
        out = (1 - min(l, dist.maxlength) * dist.p) * out
    end
    return out
end


"""
    Levenshtein()

Creates the Levenshtein distance

The Levenshtein distance is the minimum number of operations (consisting of insertions, deletions, 
substitutions of a single character) required to change one string into the other.
"""
struct Levenshtein{V <: Union{Int, Nothing}} <: Metric
   max_dist::V
end
Levenshtein() = Levenshtein(nothing)
## Source: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
# Return max_dist + 1 if distance higher than max_dist 
# to differentiate distance equal to max_dist or not, which is important for find fctions.
function (dist::Levenshtein{T})(s1, s2) where {T}
    (s1 === missing) | (s2 === missing) && return missing
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    end
    if T <: Int
        len2 - len1 > dist.max_dist && return dist.max_dist + 1
    end
    # prefix common to both strings can be ignored
    k = common_prefix(s1, s2)
    k == len1 && return len2 - k
    # distance initialized to first row of matrix
    # distance between "" and s2[1:i]
    v = collect(1:(len2-k))
    current = 0
    for (i1, ch1) in enumerate(s1)
        i1 > k || continue
        left = current = i1 - k - 1
        if T <: Int
            value_lb = left - 1
        end
        for (i2, ch2) in enumerate(s2)
            i2 > k || continue
            @inbounds above, current, left = current, left, v[i2 - k]
            if ch1 != ch2
                current = min(current, above, left) + 1
            end
            if T <: Int
                value_lb = min(value_lb, left)
            end
            @inbounds v[i2 - k] = current
        end
        if T <: Int
            value_lb > dist.max_dist && return dist.max_dist + 1
        end
    end
    if T <: Int
        current > dist.max_dist && return dist.max_dist + 1 
    end
    return current
end

"""
        OptimalStringAlignement()

    Creates the OptimalStringAlignement distance (also known ad the unrestricted DamerauLevenshtein distance).

    It is the minimum number of operations (consisting of insertions, 
    deletions or substitutions of a single character, or transposition of two adjacent characters) 
    required to change one string into the other.

    The distance differs slightly from the Damerau-Levenshtein algorithm by imposing 
    the restriction that no substring is edited more than once. So for example, "CA" to "ABC" has an edit 
    distance of 2 by a complete application of Damerau-Levenshtein, but a distance of 3 by this method that
    uses the optimal string alignment algorithm. In particular, the restricted distance does not satisfy 
    the triangle inequality.
"""
struct OptimalStringAlignement{V <: Union{Int, Nothing}} <: SemiMetric
   max_dist::V
end
OptimalStringAlignement() = OptimalStringAlignement(nothing)

## http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
# Return max_dist + 1 if distance higher than max_dist
function (dist::OptimalStringAlignement{T})(s1, s2) where {T}
    (s1 === missing) | (s2 === missing) && return missing
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    end
    if T <: Int 
        len2 - len1 > dist.max_dist && return dist.max_dist + 1
    end
    k = common_prefix(s1, s2)
    k == len1 && return len2 - k
    v = collect(1:(len2 - k))
    w = similar(v)
    prevch1, prevch2 = first(s1), first(s2)
    current = 0
    if T <: Int
        i2_start = 0
        i2_end = dist.max_dist
    end
    for (i1, ch1) in enumerate(s1)
        i1 > k || (prevch1 = ch1 ; continue)
        left = i1 - k - 1
        current = left + 1
        nextTransCost = 0
        if T <: Int
            i2_start += i1 - k - 1 + len2 - len1 > dist.max_dist
            i2_end += i2_end < len2
        end
        for (i2, ch2) in enumerate(s2)
            i2 > k || (prevch2 = ch2 ; continue)
            # no need to look beyond window of lower right diagonal - max distance cells 
            # lower right diag is i1 - (len2 - len1)) and the upper left diagonal + dist.max_dist cells (upper left is i1)
            if T <: Int
                if !(i2_start <= i2 - k - 1 < i2_end)
                    prevch2 = ch2
                    continue
                end
            end
            @inbounds above, current, left = current, left, v[i2 - k]
            @inbounds w[i2 - k], nextTransCost, thisTransCost = current, w[i2 - k], nextTransCost
            if ch1 != ch2
                current = min(current, above, left) + 1
                if i1 - k - 1 > 0 && i2 - k - 1 > 0 && ch1 == prevch2 && prevch1 == ch2
                    thisTransCost += 1
                    current = min(current, thisTransCost)
                end
            end
            @inbounds v[i2 - k] = current
            prevch2 = ch2
        end
        if T <: Int
            v[i1 - k + len2 - len1] > dist.max_dist && return dist.max_dist + 1
        end
        prevch1 = ch1
    end
    if T <: Int
        current > dist.max_dist && return dist.max_dist + 1
    end
    return current
end


"""
    DamerauLevenshtein()

Creates the DamerauLevenshtein distance

The DamerauLevenshtein distance is the minimum number of operations (consisting of insertions, 
deletions or substitutions of a single character, or transposition of two adjacent characters) 
required to change one string into the other.
"""
struct DamerauLevenshtein <: Metric end

# https://en.wikipedia.org/wiki/Damerau–Levenshtein_distance
# https://www.lemoda.net/text-fuzzy/damerau-levenshtein/
# Compared to Levenshtein & Restricted distance, cannot get by with only two vectors since transposition can be global
function (dist::DamerauLevenshtein)(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    end
    # da[ch1] will store last spotted position of ch1 in s1
    da = Dict{eltype(s1), UInt32}()
    sizehint!(da, len1)
    # distm[i1+1, i2+1] will store the distance between Iterators.take(s1, i1) and Iterators.take(s2, i2)
    distm = zeros(UInt32, len1 + 1, len2 + 1)
    distm[:, 1] = 0:len1
    distm[1, :] = 0:len2
    for (i1, ch1) in enumerate(s1)
        # j2 is last spotted position of ch1 in s2
        # j1 will be last spotted position of ch2 in s1
        j2 = 0
        for (i2, ch2) in enumerate(s2)
            match = ch1 == ch2
            @inbounds pre = min(distm[i1, i2] + !match, 
                                distm[i1 + 1, i2] + 1,
                                distm[i1, i2 + 1] + 1)
            # avoid lookup if we already know transposition won't be chosen
            j1 = (i1 == 1 || j2 == 0 || match) ? 0 : get(da, ch2, 0)
            @inbounds distm[i1 + 1, i2 + 1] = (j1 == 0) ? pre : min(pre, distm[j1, j2] + (i1 - j1 - 1) + 1 + (i2 - j2 - 1))
            if match
                j2 = i2
            end
        end
        da[ch1] = i1
    end
    return Int(distm[end, end])
end

"""
    RatcliffObershelp()

Creates the RatcliffObershelp distance

The distance between two strings is defined as one minus  the number of matching characters 
divided by the total number of characters in the two strings. Matching characters are those 
in the longest common subsequence plus, recursively, matching characters in the unmatched 
region on either side of the longest common subsequence.
"""
struct RatcliffObershelp <: SemiMetric end

function (dist::RatcliffObershelp)(s1, s2)
    (s1 === missing) | (s2 === missing) && return missing
    len1, len2 = length(s1), length(s2)
    n_matched = length_matching_blocks(s1, s2, 1, 1, len1, len2)
    len1 + len2 == 0 ? 0.0 : 1 - 2 * n_matched / (len1 + len2)
end

function length_matching_blocks(s1, s2, start1::Integer, start2::Integer, end1::Integer, end2::Integer)
    # p is just a storage vector which will be reused
    p = zeros(Int, max(end1 - start1, end2 - start2) + 1)
    length_matching_blocks!(p, s1, s2, start1, start2, end1, end2)
end

function length_matching_blocks!(p::Vector{Int}, s1, s2, start1::Integer, start2::Integer, end1::Integer, end2::Integer)
    end1 >= start1 || return 0
    end2 >= start2 || return 0
    j1, j2, len = longest_common_pattern!(p, s1, s2, start1, start2, end1, end2)
    # exit if there is no common substring
    len == 0 && return 0
    return len + 
    length_matching_blocks!(p, s1, s2, start1, start2, j1 - 1, j2 - 1) +
    length_matching_blocks!(p, s1, s2, j1 + len, j2 + len, end1, end2)
end

function longest_common_pattern!(p, s1, s2, start1, start2, end1, end2)
    if end1 - start1 > end2 - start2
        j2, j1, len = longest_common_pattern!(p, s2, s1, start2, start1, end2, end1)
    else
        j1, j2, len = 0, 0, 0
        fill!(p, 0)
        # p[i2-start2+1] stores the startingindex of the longest 
        # common pattern up to i2 with prevch1 as last matching character
        for (i1, ch1) in enumerate(s1)
            i1 >= start1 || continue
            i1 <= end1 || break
            oldj2 = 0
            for (i2, ch2) in enumerate(s2)
                i2 >= start2 || continue
                i2 <= end2 || break
                if ch1 != ch2
                    newj2 = 0
                else
                    newj2 = oldj2 > 0 ? oldj2 : i2
                    newlen = i2 - newj2 + 1
                    if newlen > len
                        j1, j2, len = i1 - newlen + 1, newj2, newlen
                    end
                end
                p[i2 - start2 + 1], oldj2 = newj2, p[i2 - start2 + 1]
            end
        end
    end
    return j1, j2, len
end
