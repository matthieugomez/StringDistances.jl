##############################################################################
##
## Normalized
##
##############################################################################

type Normalized{T}
	dist::T
end

function evaluate(normalized::Normalized, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    evaluate(normalized.dist, s1, s2, len1, len2)
end

function evaluate{T <: Union{Hamming, Levenshtein, DamerauLevenshtein}}(
	normalized::Normalized{T}, s1::AbstractString, s2::AbstractString,
    len1::Integer, len2::Integer)
    distance = evaluate(normalized.dist, s1, s2, len1, len2)
    return distance / len2
end

function evaluate{T <: QGram}(normalized::Normalized{T}, s1::AbstractString, s2::AbstractString, 
    len1::Integer, len2::Integer)
    distance = evaluate(normalized.dist, s1, s2, len1, len2)
    len1 <= (normalized.dist.q - 1)
    	return s1 == s2 ? 0.0 : 1.0
    else 
    	return distance / (len1 + len2 - 2 * normalized.dist.q + 2)
    end
end