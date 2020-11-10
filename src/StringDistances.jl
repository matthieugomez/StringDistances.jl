module StringDistances

using Distances

include("distances/utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("normalize.jl")

const StringDistance = Union{Hamming, Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp, QGramDistance, Winkler, Partial, TokenSort, TokenSet, TokenMax, Normalize}
# Distances API
Distances.result_type(dist::StringDistance, s1::Type, s2::Type) = typeof(dist("", ""))
Distances.result_type(dist::StringDistance, s1, s2) = result_type(dist, typeof(s1), typeof(s2))


include("find.jl")
include("pairwise.jl")

##############################################################################
##
## Export
##
##############################################################################

export
StringDistance,
Hamming,
Levenshtein,
DamerauLevenshtein,
Jaro,
RatcliffObershelp,
QGramDistance,
QGram,
Cosine,
Jaccard,
SorensenDice,
Overlap,
MorisitaOverlap,
NMD,
QGramDict,
QGramSortedVector,
Winkler,
Partial,
TokenSort,
TokenSet,
TokenMax,
evaluate,
compare,
result_type,
qgrams,
normalize,
findnearest,
pairwise,
pairwise!
end

