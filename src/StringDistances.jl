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
JaroWinkler, jaro_winkler, jaro,
DamerauLevenshtein, damerau_levenshtein,
QGram, qgram,
Cosine, cosine,
Jaccard, jaccard



include("edit_distances.jl")
include("qgrams_distances.jl")


end 