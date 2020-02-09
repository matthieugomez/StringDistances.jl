"""
    Jaro()

Creates the Jaro metric

The Jaro distance is defined as


``1 - (m / |s1| + m / |s2| + (m - t) / m) / 3``

where ``m`` is the number of matching characters and 
``t`` is half the number of transpositions.
"""
struct Jaro <: SemiMetric end
isnormalized(::Jaro) = true

## http://alias-i.com/lingpipe/docs/api/com/aliasi/spell/JaroWinklerDistance.html
## accepts any iterator, including AbstractString
function evaluate(dist::Jaro, s1, s2, max_dist = nothing)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    # if both are empty, m = 0 so should be 1.0 according to wikipedia. 
    # Add this line so that not the case
    len2 == 0 && return 0.0
    maxdist = max(0, div(len2, 2) - 1)
    flag = fill(false, len2)
    prevstate1 = firstindex(s1)
    i1_match = fill(prevstate1, len1)
    #  m counts number matching characters
    m = 0 
    i1 = 1
    i2 = 1
    x1 = iterate(s1)
    x2 = iterate(s2)
    while x1 !== nothing
        ch1, state1 = x1
        if i2 <= i1 - maxdist - 1
            ch2, state2 = x2
            i2 += 1
            x2 = iterate(s2, state2)
        end 
        i2curr = i2
        x2curr = x2
        while x2curr !== nothing
            i2curr > i1 + maxdist && break
            ch2, state2 = x2curr
            if (ch1 == ch2) && !flag[i2curr] 
                m += 1
                flag[i2curr] = true
                i1_match[m] = prevstate1
                break
            end
            x2curr = iterate(s2, state2) 
            i2curr += 1
        end
        x1 = iterate(s1, state1)
        i1 += 1
        prevstate1 = state1
    end
    m == 0 && return 1.0
    # t counts number of transpositions
    t = 0
    i1 = 0
    i2 = 0
    for ch2 in s2
        i2 += 1
        if flag[i2]
            i1 += 1
            t += ch2 != iterate(s1, i1_match[i1])[1]
        end
    end
    return 1.0 - (m / len1 + m / len2 + (m - t/2) / m) / 3.0
end

"""
    Levenshtein()

Creates the Levenshtein metric

The Levenshtein distance is the minimum number of operations (consisting of insertions, deletions, 
substitutions of a single character) required to change one string into the other.
"""
struct Levenshtein <: Metric end

## Source: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
# Return max_dist +1 if distance higher than max_dist
# This makes it possible to differentiate distance equalt to max_dist vs strictly higher
# This is important for find_all
## accepts any iterator, including AbstractString
function evaluate(dist::Levenshtein, s1, s2, max_dist = nothing)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    max_dist !== nothing && len2 - len1 > max_dist && return max_dist + 1
    # prefix common to both strings can be ignored
    k, x1, x2start = common_prefix(s1, s2)
    x1 == nothing && return len2 - k
    # distance initialized to first row of matrix
    # => distance between "" and s2[1:i}
    v = collect(1:(len2-k))
    current = 0
    i1 = 1
    while x1 !== nothing
        ch1, state1 = x1
        left = i1 - 1
        current = i1 - 1
        min_dist = i1 - 2 
        i2 = 1
        x2 = x2start
        while x2 !== nothing
            ch2, state2 = x2
            #  update
            above, current, left = current, left, v[i2]
            if ch1 != ch2
                current = min(current + 1, above + 1, left + 1)
            end
            min_dist = min(min_dist, left)
            v[i2] = current
            x2 = iterate(s2, state2)
            i2 += 1
        end
        max_dist !== nothing && min_dist > max_dist && return max_dist + 1
        x1 = iterate(s1, state1)
        i1 += 1
    end
    max_dist !== nothing && current > max_dist && return max_dist + 1 
    return current
end



"""
    DamerauLevenshtein()

Creates the DamerauLevenshtein metric

The DamerauLevenshtein distance is the minimum number of operations (consisting of insertions, 
deletions or substitutions of a single character, or transposition of two adjacent characters) 
required to change one string into the other.
"""
struct DamerauLevenshtein <: SemiMetric end

## http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
## accepts any iterator, including AbstractString
function evaluate(dist::DamerauLevenshtein, s1, s2, max_dist = nothing)
    (ismissing(s1) | ismissing(s2)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    max_dist !== nothing && len2 - len1 > max_dist && return max_dist + 1
    # prefix common to both strings can be ignored
    k, x1, x2start = common_prefix(s1, s2)
    (x1 == nothing) && return len2 - k
    v = collect(1:(len2-k))
    w = similar(v)
    if max_dist !== nothing
        offset = 1 + max_dist - (len2 - len1)
        i2_start = 1
        i2_end = max_dist
    end
    i1 = 1
    current = i1
    prevch1, = x1
    while x1 !== nothing
        ch1, state1 = x1
        left = (i1 - 1) 
        current = i1 
        nextTransCost = 0
        prevch2, = x2start
        if max_dist !== nothing
            i2_start += (i1 > offset) ? 1 : 0
            i2_end = min(i2_end + 1, len2)
        end
        x2 = x2start
        i2 = 1
        while x2 !== nothing
            ch2, state2 = x2
            if max_dist == nothing || (i2_start <= i2 <= i2_end)
                above = current
                thisTransCost = nextTransCost
                nextTransCost = w[i2]
                # cost of diagonal (substitution)
                w[i2] = current = left
                # left now equals current cost (which will be diagonal at next iteration)
                left = v[i2]
                if ch1 != ch2
                    # insertion
                    if left < current
                        current = left
                    end
                    # deletion
                    if above < current
                        current = above
                    end
                    current += 1
                    if (i1 != 1) & (i2 != 1) & (ch1 == prevch2) & (prevch1 == ch2)
                        thisTransCost += 1
                        if thisTransCost < current
                            current = thisTransCost
                        end
                    end
                end
                v[i2] = current
            end
            x2 = iterate(s2, state2)
            i2 += 1
            prevch2 = ch2
        end
        max_dist !== nothing && v[i1 + len2 - len1] > max_dist && return max_dist + 1
        x1 = iterate(s1, state1)
        i1 += 1
        prevch1 = ch1
    end
    max_dist !== nothing && current > max_dist && return max_dist + 1
    return current
end

"""
    RatcliffObershelp()

Creates the RatcliffObershelp metric

The distance between two strings is defined as one minus  the number of matching characters 
divided by the total number of characters in the two strings. Matching characters are those 
in the longest common subsequence plus, recursively, matching characters in the unmatched 
region on either side of the longest common subsequence.
"""
struct RatcliffObershelp <: SemiMetric end

isnormalized(::RatcliffObershelp) = true

function evaluate(dist::RatcliffObershelp, s1, s2, max_dist = nothing)
    (ismissing(s1) | ismissing(s2)) && return missing
    n_matched = sum(last.(matching_blocks(s1, s2)))
    len1, len2 = length(s1), length(s2)
    len1 + len2 == 0 ? 0. : 1.0 - 2 *  n_matched / (len1 + len2)
end

function matching_blocks(s1, s2)
    matching_blocks!(Set{Tuple{Int, Int, Int}}(), s1, s2, length(s1), length(s2), 1, 1)
end

function matching_blocks!(x::Set{Tuple{Int, Int, Int}}, s1, s2, 
    len1::Integer, len2::Integer, start1::Integer, start2::Integer)
    a = longest_common_pattern(s1, s2, len1 , len2)
    # exit if there is no common substring
    a[3] == 0 && return x
    # add the info of the common to the existing set
    push!(x, (a[1] + start1 - 1, a[2] + start2 - 1, a[3]))
    # add the longest common substring that happens before
    s1before = _take(s1, a[1] - 1)
    s2before = _take(s2, a[2] - 1)
    matching_blocks!(x, s1before, s2before, a[1] - 1, a[2] - 1, start1, start2)
    # add the longest common substring that happens after
    s1after = _drop(s1, a[1] + a[3] - 1)
    s2after = _drop(s2, a[2] + a[3] - 1)
    matching_blocks!(x, s1after, s2after, len1 - (a[1] + a[3]) + 1, 
        len2 - (a[2] + a[3]) + 1, start1 + a[1] + a[3] - 1, start2 + a[2] + a[3] - 1)
    return x
end

function longest_common_pattern(s1, s2, len1::Integer, len2::Integer)
    if len1 > len2
        start2, start1, len = longest_common_pattern(s2, s1, len2, len1)
    else
        start1, start2, len = 0, 0, 0
        p = zeros(Int, len2)
        i1 = 0
        for ch1 in s1
            i1 += 1
            i2 = 0
            oldp = 0
            for ch2 in s2
                i2 += 1
                newp = 0
                if ch1 == ch2
                    newp = oldp > 0 ? oldp : i2
                    currentlength = (i2 - newp + 1)
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