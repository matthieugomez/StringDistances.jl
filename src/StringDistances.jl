module StringDistances

using Distances: Distances, SemiMetric, Metric, evaluate, result_type
using StatsAPI: StatsAPI, pairwise, pairwise!
# Distances API
abstract type StringSemiMetric <: SemiMetric end
abstract type StringMetric <: Metric end
const StringDistance = Union{StringSemiMetric, StringMetric}
function Distances.result_type(dist::Union{StringSemiMetric, StringMetric}, s1::Type, s2::Type)
    T = typeof(dist("", ""))
    if (Missing <: s1) | (Missing <: s2)
        T = Union{T, Missing}
    end
    return T
end
Distances.result_type(dist::Union{StringSemiMetric, StringMetric}, s1, s2) = result_type(dist, typeof(s1), typeof(s2))



(dist::Union{StringSemiMetric, StringMetric})(s1, s2; max_dist = nothing) = dist(s1, s2)
include("utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("pairwise.jl")
include("normalize.jl")
include("find.jl")
include("fuzzywuzzy.jl")


##############################################################################
##
## Export
##
##############################################################################

export
StringDistance, 
StringSemiMetric,
StringMetric,
# edit distances
Hamming,
Jaro,
JaroWinkler,
Levenshtein,
OptimalStringAlignment,
DamerauLevenshtein,
RatcliffObershelp,
# Qgram distances
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
qgrams,
# normalize
compare,
# fuzzywuzzy
Partial,
TokenSort,
TokenSet,
TokenMax,
# find
findnearest,
# re-rexport from Distances.jl
evaluate,
result_type,
pairwise,
pairwise!
end

