module StringDistances

using Distances
import Distances: evaluate, result_type


abstract type StringDistance <: SemiMetric end
include("utils.jl")
include("edit.jl")
include("qgram.jl")
include("compare.jl")
include("find.jl")

##############################################################################
##
## Handle missing values
##
##############################################################################

evaluate(::QGramDistance, ::Missing, ::AbstractString) = missing
evaluate(::QGramDistance, ::AbstractString, ::Missing) = missing

evaluate(::RatcliffObershelp, ::Missing, ::AbstractString) = missing
evaluate(::RatcliffObershelp, ::AbstractString, ::Missing) = missing



function result_type(dist::StringDistance, s1::AbstractString, s2::AbstractString)
    typeof(evaluate(dist, oneunit(s1), oneunit(s2)))
end

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
qgrams
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
