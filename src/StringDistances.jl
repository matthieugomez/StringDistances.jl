__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################

import Distances: evaluate, Hamming, hamming, PreMetric, SemiMetric
import Iterators: chain
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
TokenSet

include("distances/evaluate.jl")
include("distances/edit.jl")
include("distances/qgram.jl")
include("distances/RatcliffObershelp.jl")

include("modifiers/compare.jl")
include("modifiers/winkler.jl")
include("modifiers/tokenize.jl")
include("modifiers/partial.jl")


end 