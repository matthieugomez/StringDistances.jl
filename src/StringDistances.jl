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
JaroWinkler, jaro_winkler,
QGram, qgram,
Cosine, cosine,
Jaccard, jaccard



include("edit_distances.jl")
include("qgrams_distances.jl")


end 