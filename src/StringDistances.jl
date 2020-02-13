module StringDistances

using Distances

include("utils.jl")
include("edit.jl")
include("qgram.jl")
include("normalize.jl")

const StringDistance = Union{Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp, QGramDistance, Winkler, Partial, TokenSort, TokenSet, TokenMax, Normalize}
Distances.result_type(dist::StringDistance, s1, s2) = typeof(dist("", ""))

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
compare(s1, s2, dist::StringDistance; min_score = 0.0) = 1 - normalize(dist)(s1, s2, 1 - min_score)

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

