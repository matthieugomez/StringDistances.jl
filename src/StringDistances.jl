__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################

import Distances: evaluate, Hamming, hamming
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

# Two reasons for this method
# 1. only do the switch once (becomes more complicated when calling Normalized or Winkler Adjusted distance
# 2. precomputes length(s1), length(s2) since costly for Unicode (I think it needs to traverse the array?)
function evaluate(dist, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
    	return evaluate(dist, s2, s1, len2, len1)
    else
	    return evaluate(dist, s1, s2, len1, len2)
	end
end



include("edit.jl")
include("qgram.jl")
include("normalized.jl")
include("winkler.jl")


end 