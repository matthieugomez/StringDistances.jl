module StringDistances

using Distances
import Distances: evaluate, result_type

include("utils.jl")
include("edit.jl")
include("qgram.jl")
include("modifiers.jl")

const StringDistance = Union{Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp, QGramDistance, Winkler, Partial, TokenSort, TokenSet, TokenMax}
result_type(dist::StringDistance, s1, s2) =  typeof(evaluate(dist, "", ""))

include("find.jl")

##############################################################################
##
## Export
##
##############################################################################

export
StringDistance,
Levenshtein,
DamerauLevenshtein,
Jaro,
RatcliffObershelp,
QGram,
Cosine,
Jaccard,
SorensenDice,
Overlap,
Winkler,
Partial,
TokenSort,
TokenSet,
TokenMax,
evaluate,
compare,
result_type,
qgrams,
normalize
end

