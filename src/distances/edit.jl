"""
    Hamming()

Creates the Hamming distance

The Hamming distance is defined as the number of characters that do not match
"""
struct Hamming{V <: Union{Integer, Nothing}} <: SemiMetric
   max_dist::V
end
Hamming() = Hamming(nothing)

function (dist::Hamming)(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    out = abs(length(s2) - length(s1))
    dist.max_dist !== nothing && out > dist.max_dist && return dist.max_dist + 1
    for (ch1, ch2) in zip(s1, s2)
        out += ch1 != ch2
        dist.max_dist !== nothing && out > dist.max_dist && return dist.max_dist + 1
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
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    # If both are empty, the formula in Wikipedia gives 0
    # Add this line so that not the case
    len2 == 0 && return 0.0
    d = max(0, div(len2, 2) - 1)
    flag = fill(false, len2)
    ch1_match = Vector{eltype(s1)}()
    for (i1, ch1) in enumerate(s1)
        for (i2, ch2) in enumerate(s2)
            # greedy alignement
            if (i2 <= i1 + d) && (i2 >= i1 - d) && (ch1 == ch2) && !flag[i2] 
                flag[i2] = true
                push!(ch1_match, ch1)
                break
            end
        end
    end
    #  m counts number matching characters
    m = length(ch1_match)
    m == 0 && return 1.0
    # t counts number transpositions
    t = 0
    i1 = 0
    for (i2, ch2) in enumerate(s2)
        if flag[i2]
            i1 += 1
            t += ch2 != ch1_match[i1]
        end
    end
    return 1.0 - (m / len1 + m / len2 + (m - t/2) / m) / 3.0
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
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
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
struct Levenshtein{V <: Union{Integer, Nothing}} <: Metric
   max_dist::V
end
Levenshtein() = Levenshtein(nothing)
## Source: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
# Return max_dist + 1 if distance higher than max_dist 
# to differentiate distance equal to max_dist or not, which is important for find fctions.
function (dist::Levenshtein)(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    dist.max_dist !== nothing && len2 - len1 > dist.max_dist && return dist.max_dist + 1
    # prefix common to both strings can be ignored
    k = common_prefix(s1, s2)
    k == len1 && return len2 - k
    # distance initialized to first row of matrix
    # distance between "" and s2[1:i]
    v = collect(1:(len2-k))
    current = 0
    for (i1, ch1) in enumerate(s1)
        i1 <= k && continue
        left = current = i1 - k - 1
        dist.max_dist !== nothing && (value_lb = left - 1)
        for (i2, ch2) in enumerate(s2)
            i2 <= k && continue
            above, current, left = current, left, v[i2 - k]
            if ch1 != ch2
                current = min(current, above, left) + 1
            end
            dist.max_dist !== nothing && (value_lb = min(value_lb, left))
            v[i2 - k] = current
        end
        dist.max_dist !== nothing && value_lb > dist.max_dist && return dist.max_dist + 1
    end
    dist.max_dist !== nothing && current > dist.max_dist && return dist.max_dist + 1 
    return current
end

"""
    DamerauLevenshtein()

Creates the restricted DamerauLevenshtein distance

The DamerauLevenshtein distance is the minimum number of operations (consisting of insertions, 
deletions or substitutions of a single character, or transposition of two adjacent characters) 
required to change one string into the other.

The restricted distance differs slightly from the classic Damerau-Levenshtein algorithm by imposing 
the restriction that no substring is edited more than once. So for example, "CA" to "ABC" has an edit 
distanceof 2 by a complete application of Damerau-Levenshtein, but a distance of 3 by this method that
uses the optimal string alignment algorithm. In particular, the restricted distance does not satisfy 
the triangle inequality.
"""

struct DamerauLevenshtein{V <: Union{Integer, Nothing}} <: SemiMetric
   max_dist::V
end
DamerauLevenshtein() = DamerauLevenshtein(nothing)

## http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
# Return max_dist + 1 if distance higher than max_dist
function (dist::DamerauLevenshtein)(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    dist.max_dist !== nothing && len2 - len1 > dist.max_dist && return dist.max_dist + 1
    # prefix common to both strings can be ignored
    k = common_prefix(s1, s2)
    k == len1 && return len2 - k
    v = collect(1:(len2-k))
    w = similar(v)
    if dist.max_dist !== nothing
        i2_start = 0
        i2_end = dist.max_dist
    end
    prevch1, prevch2 = first(s1), first(s2)
    current = 0
    for (i1, ch1) in enumerate(s1)
        i1 <= k && continue
        left = i1 - k - 1
        current = left + 1
        nextTransCost = 0
        if dist.max_dist !== nothing
            i2_start += (i1 - k - 1 > dist.max_dist - (len2 - len1)) ? 1 : 0
            i2_end += (i2_end < len2) ? 1 : 0
        end
        for (i2, ch2) in enumerate(s2)
            if i2 <= k 
                prevch2 = ch2
            elseif (dist.max_dist !== nothing) && ((i2 - k - 1 < i2_start) | (i2 - k - 1 >= i2_end))
                # no need to look beyond window of lower right diagonal - max distance cells 
                #lower right diag is i1 - (len2 - len1)) and the upper left diagonal + dist.max_dist cells (upper left is i1)
                prevch2 = ch2
            else
                above, current, left = current, left, v[i2 - k]
                w[i2 - k], nextTransCost, thisTransCost = current, w[i2 - k], nextTransCost
                # left now equals current cost (which will be diagonal at next iteration)
                if ch1 != ch2
                    current = min(left, current, above) + 1
                    # never happens at i2 = k + 1 because then the two previous characters were equal
                    if (i1 - k - 1 > 0) & (i2 - k - 1 > 0) && (ch1 == prevch2) && (prevch1 == ch2)
                        thisTransCost += 1
                        current = min(current, thisTransCost)
                    end
                end
                v[i2 - k] = current
                prevch2 = ch2
            end
        end
        dist.max_dist !== nothing && v[i1 - k + len2 - len1] > dist.max_dist && return dist.max_dist + 1
        prevch1 = ch1
    end
    dist.max_dist !== nothing && current > dist.max_dist && return dist.max_dist + 1
    return current
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
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    n_matched = sum(last.(matching_blocks(s1, s2)))
    len1, len2 = length(s1), length(s2)
    len1 + len2 == 0 ? 0. : 1.0 - 2 *  n_matched / (len1 + len2)
end

function matching_blocks(s1, s2)
    matching_blocks!(Set{Tuple{Int, Int, Int}}(), s1, s2, 1, 1)
end

function matching_blocks!(x::Set{Tuple{Int, Int, Int}}, s1, s2, start1::Integer, start2::Integer)
    n1, n2, len = longest_common_pattern(s1, s2)
    # exit if there is no common substring
    len == 0 && return x
    # add the info of the common to the existing set
    push!(x, (n1 + start1 - 1, n2 + start2 - 1, len))
    # add the longest common substring that happens before
    matching_blocks!(x, _take(s1, n1 - 1), _take(s2, n2 - 1), start1, start2)
    # add the longest common substring that happens after
    matching_blocks!(x, _drop(s1, n1 + len - 1), _drop(s2, n2 + len - 1), 
                    start1 + n1 + len - 1, start2 + n2 + len - 1)
    return x
end

function longest_common_pattern(s1, s2)
    if length(s1) > length(s2)
        start2, start1, len = longest_common_pattern(s2, s1)
    else
        start1, start2, len = 0, 0, 0
        p = zeros(Int, length(s2))
        for (i1, ch1) in enumerate(s1)
            oldp = 0
            for (i2, ch2) in enumerate(s2)
                newp = 0
                if ch1 == ch2
                    newp = oldp > 0 ? oldp : i2
                    currentlength = i2 - newp + 1
                    if currentlength > len
                        start1, start2, len = i1 - currentlength + 1, newp, currentlength
                    end
                end
                p[i2], oldp = newp, p[i2]
            end
        end
    end
    return start1, start2, len
end
