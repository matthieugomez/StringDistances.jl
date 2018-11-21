
module StringDistances

##############################################################################
##
## Export
##
##############################################################################
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

