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
## TypeAlias
##
##############################################################################

const GraphemeIterator = Base.UTF8proc.GraphemeIterator
#const AbstractStringorGraphemeIterator = Union{AbstractString, Base.UTF8proc.GraphemeIterator}
const AbstractStringorGraphemeIterator = AbstractString

##############################################################################
##
## GraphemeIterator. Unicode iteration broken because Unicode 9 has different iteration properties.
##
##############################################################################
#Base.prevind(x::GraphemeIterator, i::Integer) = prevind(x.s, i)
#Base.nextind(x::GraphemeIterator, i::Integer) = nextind(x.s, i)
#Base.chr2ind(x::GraphemeIterator, i::Integer) = chr2ind(x.s, i)
#Base.SubString(x::GraphemeIterator, i::Integer, j::Integer) = graphemeiterator(SubString(x.s, i::Integer, j:#:Integer))
#graphemeiterator(s::AbstractString) = GraphemeIterator{typeof(s)}(s)
#
## added
##these 2 functions allow to define prevind nextind, chr2ind, prevind etc
#function Base.isvalid(s::GraphemeIterator, i::Integer)
#    if !isvalid(s.s, i)
#        return false
#    else
#        k = start(s)
#        while !done(s, k)
#            j = k[1]
#            if j == i
#                return true
#            end
#        end
#        return false
#    end
#end
#function Base.endof(s::GraphemeIterator)
#    k = start(s)
#    while !done(s, k)
#        i = k[1]
#        c, k = next(s, k)
#    end
#    return i
#end
#
## 1. issues with stuff like search  or print_escaped where character iteration vs string iteration matters. #I need to pass the original string for now
#Base.search(x::GraphemeIterator, s::Vector{Char}) = search(x.s, s)
## 2. issue with keeping iterator property for stuff like split, join. for now, I decide to loose the #enumerator property but add it back after join. But SubString for instance does not loose the property
#Base.split(x::GraphemeIterator, args...) = split(x.s, args...)
#iterator{T <: GraphemeIterator}(::Type{T}, x::AbstractString) = graphemeiterator(x)
iterator(::Type{T}, x::AbstractString) where {T <: AbstractString} = x
#
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
for x in (:evaluate, :compare)
    @eval begin
        function $x(dist::PreMetric, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator)
            len1, len2 = length(s1), length(s2)
            if len1 > len2
                return $x(dist, s2, s1, len2, len1)
            else
                return $x(dist, s1, s2, len1, len2)
            end
        end
    end
end


##############################################################################
##
## compare
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator,
    len1::Integer, len2::Integer)
    1.0 - evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein},
    s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator,
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    len2 == 0 ? 1.0 : 1.0 - distance / len2
end

# compare always return a value between 0 and 1.
# When string length < q for qgram distance, returns s1 == s2
function compare(dist::AbstractQGram,
    s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator,
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::QGram,
    s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator,
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    distance = evaluate(dist, s1, s2, len1, len2)
    1 - distance / (len1 + len2 - 2 * dist.q + 2)
end





end
