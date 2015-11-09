__precompile__(true)

module StringDistances

##############################################################################
##
## Export
##
##############################################################################
import Base.UTF8proc.GraphemeIterator
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
##############################################################################
##
## Extend methods of GraphemeIterator
##
##############################################################################
typealias GraphemeOrString Union{GraphemeIterator, AbstractString}

# add the following methods
Base.nextind(g::GraphemeIterator, state::Integer) = next(g, state)[2]
function Base.chr2ind(s::GraphemeIterator, i::Integer)
    i < start(s) && throw(BoundsError(s.s, i))
    j = 1
    k = start(s)
    while true
        c, l = next(s,k)
        if i == j
            return k
        end
        j += 1
        k = l
    end
end
Base.endof(g::GraphemeIterator) = endof(g.s)
Base.SubString(x::GraphemeIterator, i, j) = graphemes(SubString(x.s, i, j))

##############################################################################
##
## include
##
##############################################################################
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

function evaluate(dist::PreMetric, s1::GraphemeOrString, s2::GraphemeOrString)
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

function compare(dist::PreMetric, s1::GraphemeOrString, s2::GraphemeOrString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        return compare(dist, s2, s1, len2, len1)
    else
        return compare(dist, s1, s2, len1, len2)
    end
end

function compare(dist::PreMetric, s1::GraphemeOrString, s2::GraphemeOrString, 
    len1::Integer, len2::Integer)
    1.0 - evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::GraphemeOrString, s2::GraphemeOrString,
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    len2 == 0 ? 1.0 : 1.0 - distance / len2
end

# compare always return a value between 0 and 1. 
# When string length < q for qgram distance, returns s1 == s2
function compare(dist::AbstractQGram, 
    s1::GraphemeOrString, s2::GraphemeOrString, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::QGram, 
    s1::GraphemeOrString, s2::GraphemeOrString, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    distance = evaluate(dist, s1, s2, len1, len2)
    1 - distance / (len1 + len2 - 2 * dist.q + 2)
end





end 