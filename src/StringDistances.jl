__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################

import Distances: evaluate, Hamming, hamming, PreMetric, SemiMetric
export evaluate,
Hamming, hamming,
Levenshtein, levenshtein,
DamerauLevenshtein, damerau_levenshtein,
Jaro, jaro,
QGram, qgram,
Cosine, cosine,
Jaccard, jaccard,
Normalized,
Winkler


# 1. only do the switch once
# 2. precomputes length(s1), length(s2)
function evaluate(dist::PreMetric, s1::AbstractString, s2::AbstractString, x...)
	len1, len2 = length(s1), length(s2)
	if len1 > len2
		return evaluate(dist, s2, s1, len2, len1, x...)
	else
		return evaluate(dist, s1, s2, len1, len2, x...)
	end
end

include("edit.jl")
include("qgram.jl")
include("normalized.jl")
include("winkler.jl")


end 