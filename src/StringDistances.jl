__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################

import Distances: evaluate

export Hamming,
Levenshtein,
JaroWinkler,
DamerauLevenshtein,
hamming,
levenshtein,
damerau_levenshtein,
jaro_winkler,
jaro

##############################################################################
##
## Hamming
##
##############################################################################
type Hamming end

function evaluate(dist::Hamming, s1::AbstractString, s2::AbstractString)
	length(s1) > length(s2) && return evaluate(dist, s2, s1)
	count = 0
	@inbounds for i in 1:length(s1)
	   count += s1[i] != s2[i]
	end
	count += length(s2) - length(s1)
	return count
end
hamming(s1::AbstractString, s2::AbstractString) = evaluate(Hamming(), s1, s2)


##############################################################################
##
## Levenshtein and Damerau Levenshtein
## Source Levenshtein: http://blog.softwx.net/2014/12/optimizing-levenshtein-algorithm-in-c.html
##  Source DamerauLevenshtein: http://blog.softwx.net/2015/01/optimizing-damerau-levenshtein_15.html
##
##############################################################################

function common_suffix(s1::AbstractString, s2::AbstractString)
    len1 = length(s1)
    len2 = length(s2)
    while ((len1 > 0) && (s1[len1] == s2[len2]))
        len1 -= 1
        len2 -= 1
    end
    return len1, len2
end

function common_prefix(s1::AbstractString, s2::AbstractString, len1::Int, len2::Int)
    start = 0
    len1 == 0 && return len1, len2, start
    if (s1[start + 1] == s2[start + 1]) 
        while ((start < len1) && (s1[start + 1] == s2[start + 1]))
            start += 1
        end
        len1 -= start
        len2 -= start
        len1 == 0 && return len1, len2, start
    end
    return len1, len2, start
end

type Levenshtein end

function evaluate(dist::Levenshtein, s1::AbstractString, s2::AbstractString)
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0

    # common 
    len1, len2 = common_suffix(s1, s2)
    len1, len2, start = common_prefix(s1, s2, len1, len2)
    len1 == 0 && return len2

    dist = Array(Int, len2)
    @inbounds for i2 in 1:len2
        dist[i2] = i2
    end
    current = 0
    for i1 in 1:len1
        ch1 = s1[start + i1]
        left = current = i1 - 1
        for i2 in 1:len2
            above = current
            current = left
            left = dist[i2]
            if ch1 != s2[start + i2]
                current += 1
                insDel = above + 1
                if insDel < current
                    current = insDel
                end
                insDel = left + 1
                if insDel < current
                    current = insDel
                end
            end
            dist[i2] = current
        end
    end
    return current
end
levenshtein(s1::AbstractString, s2::AbstractString) = evaluate(Levenshtein(), s1, s2)

type DamerauLevenshtein end

function evaluate(dist::DamerauLevenshtein, s1::AbstractString, s2::AbstractString)
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0

    # common 
    len1, len2 = common_suffix(s1, s2)
    len1, len2, start = common_prefix(s1, s2, len1, len2)
    len1 == 0 && return len2

    dist = Array(Int, length(s2))
    @inbounds for i2 in 1:len2
        dist[i2] = i2
    end
    dist2 = Array(Int, length(s2))

    ch1 = s1[1]
    current = 0
    for i1 in 1:len1
        prevch1 = ch1
        ch1 = s1[start + i1]
        ch2 = s2[start + 1]
        left = i1 - 1
        current = i1
        nextTransCost = 0
        for i2 in 1:len2
            above = current
            thisTransCost = nextTransCost
            nextTransCost = dist2[i2]
            dist2[i2] = current = left
            left = dist[i2]
            prevch2 = ch2
            ch2 = s2[start + i2]
            if ch1 != ch2
                if left < current
                    current = left
                end
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
            dist[i2] = current
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

function evaluate(dist::JaroWinkler, s1::AbstractString, s2::AbstractString) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0.0

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
    m == 0.0 && return 0.0
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
    return score
end

function jaro_winkler(s1::AbstractString, s2::AbstractString; 
        scaling_factor = 0.1, boosting_threshold = 0.7, long_threshold = 5)
    evaluate(JaroWinkler(scaling_factor, boosting_threshold, long_threshold), s1, s2)
end

jaro(s1::AbstractString, s2::AbstractString) = evaluate(JaroWinkler(0.0, 0.0, 0), s1, s2)

end 