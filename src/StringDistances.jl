module StringDistances

using Distances
import StatsAPI: pairwise, pairwise!
abstract type StringSemiMetric <: SemiMetric end
abstract type StringMetric <: Metric end
const StringDistance = Union{StringSemiMetric, StringMetric}
function Distances.result_type(dist::StringDistance, s1::Type, s2::Type)
    T = typeof(dist("", ""))
    if (Missing <: s1) | (Missing <: s2)
        T = Union{T, Missing}
    end
    return T
end
Distances.result_type(dist::StringDistance, s1, s2) = result_type(dist, typeof(s1), typeof(s2))


include("distances/utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")


include("normalize.jl")
include("pairwise.jl")
include("find.jl")
include("fuzzywuzzy.jl")



##############################################################################
##
## Export
##
##############################################################################

export StringDistance,
StringSemiMetric,
StringMetric,
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

