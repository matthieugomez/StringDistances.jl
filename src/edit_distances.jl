
##############################################################################
##
## Hamming
##
##############################################################################

function evaluate(dist::Hamming, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    len1 > len2 && return evaluate(dist, s2, s1)
    count = 0

    state2 = start(s2)
    for ch1 in s1
        ch2, state2 = next(s2, state2)
        count += ch1 != ch2
    end
    count += len2 - len1
    return count
end

hamming(s1::AbstractString, s2::AbstractString) = evaluate(Hamming(), s1, s2)

##############################################################################
##
## Levenshtein and Damerau Levenshtein
## Source Levenshtein: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
## Source DamerauLevenshtein: http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
##
##############################################################################
# prefix common to both strings can be ignored
function common_prefix(s1::AbstractString, s2::AbstractString)
    start1 = start(s1)
    start2 = start(s2)
    while !done(s1, start1)
        ch1, nextstart1 = next(s1, start1)
        ch2, nextstart2 = next(s2, start2)
        ch1 != ch2 && break
        start1, start2 = nextstart1, nextstart2
    end
    return start1, start2
end
type Levenshtein end
function evaluate(dist::Levenshtein, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)

    len1 > len2 && return evaluate(dist, s2, s1)
    len2 == 0 && return 0

    # common 
    start1, start2 = common_prefix(s1, s2)
    done(s1, start1) && return len2

    # distance initialized to first row of matrix
    # => distance between "" and s2[1:i}
    v0 = Array(Int, len2)
    @inbounds for i2 in 1:len2
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
function levenshtein(s1::AbstractString, s2::AbstractString)
    evaluate(Levenshtein(), s1, s2)
end

type DamerauLevenshtein end

function evaluate(dist::DamerauLevenshtein, s1::AbstractString, s2::AbstractString)
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0

    # common 
    len1, len2 = length(s1), length(s2)
    start1, start2 = common_prefix(s1, s2)
    done(s1, start1) && return len2

    v0 = Array(Int, length(s2))
    @inbounds for i2 in 1:len2
        v0[i2] = i2
    end
    v2 = Array(Int, length(s2))

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

damerau_levenshtein(s1::AbstractString, s2::AbstractString) = evaluate(DamerauLevenshtein(), s1, s2)

##############################################################################
##
## JaroWinkler
##
##############################################################################

type JaroWinkler{T1 <: Number, T2 <: Number, T3 <: Integer}
    scaling_factor::T1      # scaling factor. Default to 0.1
    boosting_threshold::T2      # boost threshold. Default to 0.7
    long_threshold::T3  # long string adjustment. Default to 5
end
JaroWinkler() = JaroWinkler(0.1, 0.7, 5)

function evaluate(dist::JaroWinkler, s1::AbstractString, s2::AbstractString) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 1.0

    maxdist = max(0, div(length(s2), 2) - 1)
    m = 0 # matching characters
    t = 0 # half number of transpositions
    flag = fill(false, length(s2))
    prevpos = 0
    @inbounds for i1 in 1:length(s1)
        ch = s1[i1]
        i2low =  max(1, i1 - maxdist)
        i2high = min(length(s2), i1 + maxdist)
        for i2 in i2low:i2high
            if ch == s2[i2] && !flag[i2] 
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
    score = (m / length(s1) + m / length(s2) + (m - t) / m) / 3.0

    # common prefix adjustment
    if (dist.scaling_factor > 0  && score >= dist.boosting_threshold) || (length(s1) >= dist.long_threshold)
        l = 0
        last = min(4, length(s1))
        while l < last && s1[l+1] == s2[l+1]
            l += 1
        end
        # common prefix adjustment
        if (dist.scaling_factor > 0  && score >= dist.boosting_threshold)
            score += l * (1 - score) * dist.scaling_factor
        end
        # longer string adjustment
        if (length(s1) >= dist.long_threshold) &&  (m - l >= 2) && ((m - l) >= (length(s1) - l) / 2)
            score += (1 - score) * (m - (l + 1)) / (length(s1) + length(s2) - (2 * (l - 1)))
        end
    end
    return 1 - score
end

function jaro_winkler(s1::AbstractString, s2::AbstractString; 
        scaling_factor = 0.1, boosting_threshold = 0.7, long_threshold = 5)
    evaluate(JaroWinkler(scaling_factor, boosting_threshold, long_threshold), s1, s2)
end

jaro(s1::AbstractString, s2::AbstractString) = evaluate(JaroWinkler(0.0, 0.0, 0), s1, s2)















