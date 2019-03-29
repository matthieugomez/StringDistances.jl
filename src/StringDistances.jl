module StringDistances

##############################################################################
##
## Export
##
##############################################################################


# some memo
# length: number of characters
# ncodeunits: Return the number of code units in a string (aking to index of vector). Not all such indices  are valid – they may not be the start of a character,.
# sizeof:  Size, in bytes, of the string str. Equal to the number of code units in str  multiplied by the size, in bytes, of one code unit in str.

import Base: eltype, length, iterate, ==, hash, isless, convert, show
import Distances: evaluate, Hamming, hamming, PreMetric, SemiMetric
import IterTools: chain
export
evaluate,
compare,
Hamming,
Levenshtein,
DamerauLevenshtein,
Jaro,
QGram,
Cosine,
Jaccard,
SorensenDice,
Overlap,
RatcliffObershelp,
Winkler,
Partial,
TokenSort,
TokenSet,
TokenMax

##############################################################################
##
## include
##
##############################################################################
include("utils.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("distances/RatcliffObershelp.jl")
include("compare.jl")

end

