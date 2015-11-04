##############################################################################
##
## compare
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        return compare(dist, s2, s1, len2, len1)
    else
        return compare(dist, s1, s2, len1, len2)
    end
end



function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    1.0 - evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, s1::AbstractString, s2::AbstractString,
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    return len2 == 0 ? 1.0 : 1.0 - distance / len2
end

function compare(dist::QGram, s1::AbstractString, s2::AbstractString, 
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    if len1 <= (dist.q - 1)
        return s1 == s2 ? 1.0 : 0.0
    else 
        return 1 - distance / (len1 + len2 - 2 * dist.q + 2)
    end
end