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
TokenSet,
TokenMax

include("distances/edit.jl")
include("distances/qgram.jl")
include("distances/RatcliffObershelp.jl")

include("modifiers/winkler.jl")
include("modifiers/fuzzywuzzy.jl")

##############################################################################
##
## Higher level functions
##
##############################################################################

function evaluate(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        return evaluate(dist, s2, s1, len2, len1)
    else
        return evaluate(dist, s1, s2, len1, len2)
    end
end

##############################################################################
##
## compare
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        return compare(dist, s2, s1, len2, len1)
    else
        return compare(dist, s1, s2, len1, len2)
    end
end

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString, 
    len1::Integer, len2::Integer)
    1.0 - evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::AbstractString, s2::AbstractString,
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    len2 == 0 ? 1.0 : 1.0 - distance / len2
end

# compare always return a value between 0 and 1. 
# When string length < q for qgram distance, returns s1 == s2
function compare(dist::AbstractQGram, 
    s1::AbstractString, s2::AbstractString, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::QGram, 
    s1::AbstractString, s2::AbstractString, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    distance = evaluate(dist, s1, s2, len1, len2)
    1 - distance / (len1 + len2 - 2 * dist.q + 2)
end


end 