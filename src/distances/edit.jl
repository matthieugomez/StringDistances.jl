##############################################################################
##
## Find common prefixes (up to lim. -1 means Inf)
##############################################################################

function common_prefix(s1::AbstractString, s2::AbstractString, lim::Integer = -1)
    x1 = iterate(s1)
    x2 = iterate(s2)
    l = 0
    while (x1 !== nothing) && (x2 !== nothing) && (l < lim || lim < 0)
        ch1, state1 = x1
        ch2, state2 = x2
        ch1 != ch2 && break
        x1 = iterate(s1, state1)
        x2 = iterate(s2, state2)
        l += 1
    end
    return l, x1, x2
end

##############################################################################
##
## Hamming
##
##############################################################################

function evaluate(dist::Hamming, s1::AbstractString, s2::AbstractString)
    out = 0
    for (ch1, ch2) in zip(s1, s2)
        out += ch1 != ch2
    end
    out += abs(length(s2) - length(s1))
    return out
end

##############################################################################
##
## Levenshtein
## Source: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
##
##############################################################################


struct Levenshtein <: SemiMetric end

function evaluate(dist::Levenshtein, s1::AbstractString, s2::AbstractString)
    # prefix common to both strings can be ignored
    s2, len2, s1, len1 = reorder(s1, s2)
    k, x1, x2start = common_prefix(s1, s2)
    (x1 == nothing) && return len2 - k
    # distance initialized to first row of matrix
    # => distance between "" and s2[1:i}
    v0 = collect(1:(len2 - k))
    current = 0
    i1 = 1
    while x1 !== nothing
        ch1, state1 = x1
        left = (i1 - 1)
        current = (i1 - 1)
        i2 = 1
        x2 = x2start
        while x2 !== nothing
            ch2, state2 = x2
            #  update
            above, current, left = current, left, v0[i2]
            if ch1 != ch2
                # substitution
                current = min(current + 1,
                                above + 1,
                                left + 1)
            end
            v0[i2] = current
            x2 = iterate(s2, state2)
            i2 += 1
        end
        x1 = iterate(s1, state1)
        i1 += 1
    end
    return current
end

##############################################################################
##
## Damerau Levenshtein
## Source: http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
##
##############################################################################

struct DamerauLevenshtein <: SemiMetric end

function evaluate(dist::DamerauLevenshtein, s1::AbstractString, s2::AbstractString)
    s2, len2, s1, len1 = reorder(s1, s2)
    # prefix common to both strings can be ignored
    k, x1, x2start = common_prefix(s1, s2)
    (x1 == nothing) && return len2 - k
    v0 = collect(1:(len2 - k))
    v2 = similar(v0)
    i1 = 1
    current = i1
    prevch1, = x1
    while (x1 !== nothing)
        ch1, state1 = x1
        left = (i1 - 1) 
        current = i1 
        nextTransCost = 0
        prevch2, = x2start
        x2 = x2start
        i2 = 1
        while (x2 !== nothing)
            ch2, state2 = x2
            above = current
            thisTransCost = nextTransCost
            nextTransCost = v2[i2]
            # cost of diagonal (substitution)
            v2[i2] = current = left
            # left now equals current cost (which will be diagonal at next iteration)
            left = v0[i2]
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
                if i1 != 1 && i2 != 1 && ch1 == prevch2 && prevch1 == ch2
                    thisTransCost += 1
                    if thisTransCost < current
                        current = thisTransCost
                    end
                end
            end
            v0[i2] = current
            x2 = iterate(s2, state2)
            i2 += 1
            prevch2 = ch2
        end
        x1 = iterate(s1, state1)
        i1 += 1
        prevch1 = ch1
    end
    return current
end

##############################################################################
##
## Jaro
## http://alias-i.com/lingpipe/docs/api/com/aliasi/spell/JaroWinklerDistance.html
##############################################################################

struct Jaro <: SemiMetric end

function evaluate(dist::Jaro, s1::AbstractString, s2::AbstractString)
    s2, len2, s1, len1 = reorder(s1, s2)
    # if both are empty, m = 0 so should be 1.0 according to wikipedia. Add this line so that not the case
    len2 == 0 && return 0.0
    maxdist = max(0, div(len2, 2) - 1)
    flag = fill(false, len2)
    prevstate1 = firstindex(s1)
    i1_match = prevstate1 * ones(Int, len1)
    #  m counts matching characters
    m = 0 
    i1 = 1
    i2 = 1
    x1 = iterate(s1)
    x2 = iterate(s2)
    while (x1 !== nothing)
        ch1, state1 = x1
        if i2 <= i1 - maxdist - 1
            ch2, state2 = x2
            i2 += 1
            x2 = iterate(s2, state2)
        end 
        i2curr = i2
        x2curr = x2
        while (x2curr !== nothing)
            (i2curr > i1 + maxdist) && break
            ch2, state2 = x2curr
            if ch1 == ch2 && !flag[i2curr] 
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
    # count t transpotsitions
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
    m == 0 && return 1.0
    score = (m / len1 + m / len2 + (m - t/2) / m) / 3.0
    return 1.0 - score
end