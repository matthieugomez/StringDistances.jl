module StringDistances

using Distances

include("utils.jl")
include("edit.jl")
include("qgram.jl")
include("modifiers.jl")

const StringDistance = Union{Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp, QGramDistance, Winkler, Partial, TokenSort, TokenSet, TokenMax, Normalize}
Distances.result_type(dist::StringDistance, s1, s2) =  typeof(dist("", ""))
Distances.evaluate(dist::StringDistance, args...) = dist(args...)
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

