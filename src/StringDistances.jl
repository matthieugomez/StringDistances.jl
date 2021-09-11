module StringDistances

using Distances
import StatsAPI: pairwise, pairwise!

include("distances/utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("normalize.jl")
include("fuzzywuzzy.jl")
const StringDistance = Union{Hamming, Jaro, JaroWinkler, Levenshtein, OptimalStringAlignement, DamerauLevenshtein, RatcliffObershelp, AbstractQGramDistance, Normalized, Partial, TokenSort, TokenSet, TokenMax}
"""
    compare(s1, s2, dist)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the distance `dist`.

### Examples
```julia-repl
julia> compare("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function compare(s1, s2, dist::StringDistance; min_score = 0.0)
    1 - normalize(dist, max_dist = 1 - min_score)(s1, s2)
end 
include("find.jl")
include("pairwise.jl")

# Distances API
function Distances.result_type(dist::StringDistance, s1::Type, s2::Type)
    T = typeof(dist("", ""))
    if (Missing <: s1) | (Missing <: s2)
        T = Union{T, Missing}
    end
    return T
end
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

