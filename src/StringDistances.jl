module StringDistances

using Distances

include("distances/utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("modifiers.jl")
include("normalize.jl")
include("pairwise.jl")
include("find_partial.jl")

# Distances API
Distances.result_type(dist::StringDistance, s1::Type, s2::Type) = typeof(dist("", ""))
Distances.result_type(dist::StringDistance, s1, s2) = result_type(dist, typeof(s1), typeof(s2))

# ambiguity fix
Distances.result_type(dist::StringDistance, s1::AbstractArray, s2::AbstractArray) = result_type(dist, eltype(s1), eltype(s2))





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
JaroWinkler,
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
findnearest_partial,
findall_partial,
pairwise,
pairwise!
end
