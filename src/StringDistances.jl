__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################
import Base: eltype, length, start, done, next, ==, hash, isless, convert, show, endof
import Base.UTF8proc: isgraphemebreak
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
longest_common_substring,
matching_blocks,
RatcliffObershelp,
Winkler,
Partial,
TokenSort,
TokenSet,
TokenMax,
graphemeiterator

##############################################################################
##
## include
##
##############################################################################
include("utils.jl")

include("distances/edit.jl")
include("distances/qgram.jl")
include("distances/RatcliffObershelp.jl")

include("modifiers/winkler.jl")
include("modifiers/fuzzywuzzy.jl")

include("compare.jl")

end

