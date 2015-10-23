__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################

import Distances: evaluate, Hamming, hamming

export Hamming,
Levenshtein,
JaroWinkler,
DamerauLevenshtein,
QGram,
Cosine,
Jaccard,
hamming,
levenshtein,
damerau_levenshtein,
jaro_winkler,
jaro,
qgram,
cosine,
jaccard


include("edit_distances.jl")
include("qgrams_distances.jl")


end 