##############################################################################
##
## Find common prefixes (up to lim. -1 means Inf)
## Assumes length(s1) <= length(s2)
##############################################################################

function common_prefix(s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, lim::Integer = -1)
    start1 = start(s1)
    start2 = start(s2)
    l = 0
    while !done(s1, start1) && (l < lim || lim < 0)
        ch1, nextstart1 = next(s1, start1)
        ch2, nextstart2 = next(s2, start2)
        ch1 != ch2 && break
        l += 1
        start1, start2 = nextstart1, nextstart2
    end
    return l, start1, start2
end

##############################################################################
##
## Hamming
##
##############################################################################

function evaluate(dist::Hamming, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, len1::Integer, len2:: Integer)
    count = 0
    for (ch1, ch2) in zip(s1, s2)
        count += ch1 != ch2
    end
    count += len2 - len1
    return count
end

##############################################################################
##
## Levenshtein
## Source: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
##
##############################################################################


type Levenshtein <: SemiMetric end
function evaluate(dist::Levenshtein, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, len1::Integer, len2::Integer)

    # prefix common to both strings can be ignored
    k, start1, start2 = common_prefix(s1, s2)
    done(s1, start1) && return len2 - k

    # distance initialized to first row of matrix
    # => distance between "" and s2[1:i}
    v0 = Array(Int, len2 - k)
    @inbounds for i2 in 1:(len2 - k)
        v0[i2] = i2 
    end
    current = zero(0)
    state1 = start1
    i1 = 0
    while !done(s1, state1)
        i1 += 1
        ch1, state1 = next(s1, state1)
        left = (i1 - 1)
        current = (i1 - 1)
        state2 = start2
        i2 = 0
        while !done(s2, state2)
            i2 += 1
            ch2, state2 = next(s2, state2)
            #  update
            above, current, left = current, left, v0[i2]
            if ch1 != ch2
                # substitution
                current = min(current + 1,
                                above + 1,
                                left + 1)
            end
            v0[i2] = current
        end
    end
    return current
end

##############################################################################
##
## Damerau Levenshtein
## Source: http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
##
##############################################################################

type DamerauLevenshtein <: SemiMetric end

function evaluate(dist::DamerauLevenshtein, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, len1::Integer, len2::Integer)

    # prefix common to both strings can be ignored
    k, start1, start2 = common_prefix(s1, s2)
    done(s1, start1) && return len2 - k

    v0 = Array(Int, len2 - k)
    @inbounds for i2 in 1:(len2 - k)
        v0[i2] = i2
    end
    v2 = Array(Int, len2 - k)

    ch1, = next(s1, start1)
    current = 0
    state1 = start1
    i1 = 0
    while !done(s1, state1)
        i1 += 1
        prevch1 = ch1
        ch1, state1 = next(s1, state1)
        ch2, = next(s2, start2)
        left = (i1 - 1) 
        current = i1 
        nextTransCost = 0
        state2 = start2
        i2 = 0
        while !done(s2, state2)
            i2 += 1
            prevch2 = ch2
            ch2, state2 = next(s2, state2)
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
        end
    end
    return current
end

##############################################################################
##
## Jaro
##
##############################################################################

type Jaro <: SemiMetric end

function evaluate(dist::Jaro, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, len1::Integer, len2::Integer) 
    # if len2 == 0, m = 0 so should be 1.0 according to wikipedia. Nope.
    len2 == 0 && return 0.0

    maxdist = max(0, div(len2, 2) - 1)
    m = 0 # matching characters
    t = 0 # half number of transpositions
    flag = fill(false, len2)
    prevpos = 0

    i1 = 0
    startstate2 = start(s2)
    starti2 = 0
    for ch1 in s1
        i1 += 1
        if starti2 < i1 - maxdist - 1
            startstate2 = nextind(s2, startstate2)
            starti2 += 1
        end 
        i2 = starti2
        state2 = startstate2
        while !done(s2, state2) && i2 < i1 + maxdist
            ch2, state2 = next(s2, state2)
            i2 += 1
            if ch1 == ch2 && !flag[i2] 
                m += 1
                # if match is before the index of previous match
                if i2 < prevpos
                    t += 1
                end
                prevpos = max(i2, prevpos)
                flag[i2] = true
                break
            end
        end
    end
    m == 0.0 && return 1.0
    score = (m / len1 + m / len2 + (m - t) / m) / 3.0
    return 1.0 - score
end

jaro(s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator) = evaluate(Jaro(), s1, s2)
