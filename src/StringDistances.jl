module StringDistances

using Distances
import Distances: evaluate, result_type

isnormalized(dist::SemiMetric) = false


include("utils.jl")
include("edit.jl")
include("qgram.jl")
include("modifier.jl")

const StringDistance = Union{Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp, QGramDistance, Winkler, Partial, TokenSort, TokenSet, TokenMax}
function result_type(dist::StringDistance, s1, s2)
    typeof(evaluate(dist, "", ""))
end

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
compare(s1, s2, dist::StringDistance; min_score = 0.0) = 1 - evaluate(normalize(dist), s1, s2, 1 - min_score)

const metrics = [Cosine, DamerauLevenshtein, Jaccard, Jaro, Levenshtein, Overlap, Partial, QGram,
    RatcliffObershelp, SorensenDice, StringDistance,  TokenMax,
    TokenSet, TokenSort, Winkler]

for M in metrics
    @eval @inline (dist::$M)(a::AbstractString, b::AbstractString) = evaluate(dist, a, b)
end

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

##############################################################################
##
## Some things about Strings

# length: number of characters
# ncodeunits: Return the number of code units in a string (aking to index of vector). 
# Not all such indices  are valid – they may not be the start of a character,.
# sizeof:  Size, in bytes, of the string str. Equal to the number of code units in str  
# multiplied by the size, in bytes, of one code unit in str.

# lastindex: Return the last index of a collection
# nextinds(s, i):  return the index of the start of the character whose encoding starts after index i
# nextind(s, 0, N): return the index of the Nth character of s (or, if there are 
# less than N characters, return ncodeunits(str) + (N - length(s))

##############################################################################
