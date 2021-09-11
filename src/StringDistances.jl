module StringDistances

using Distances
import StatsAPI: pairwise, pairwise!

include("distances/utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("normalize.jl")
include("fuzzywuzzy.jl")
const StringDistance = Union{Hamming, Jaro, JaroWinkler, Levenshtein, OptimalStringAlignement, DamerauLevenshtein, RatcliffObershelp, AbstractQGramDistance, Normalized, Partial, TokenSort, TokenSet, TokenMax}
include("compare.jl")
include("pairwise.jl")

# Distances API
Distances.result_type(dist::StringDistance, s1::Type, s2::Type) = typeof(dist("", ""))
Distances.result_type(dist::StringDistance, s1, s2) = result_type(dist, typeof(s1), typeof(s2))





##############################################################################
##
## Export
##
##############################################################################

export
StringDistance,
Hamming,
Jaro,
JaroWinkler,
Levenshtein,
OptimalStringAlignement,
DamerauLevenshtein,
RatcliffObershelp,
AbstractQGramDistance,
QGramDict,
QGramSortedVector,
QGram,
Cosine,
Jaccard,
SorensenDice,
Overlap,
MorisitaOverlap,
NMD,
Partial,
TokenSort,
TokenSet,
TokenMax,
evaluate,
compare,
result_type,
qgrams,
findnearest,
pairwise,
pairwise!
end

