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
## Define 2 iterators. One on character one on Grapheme.
##
##############################################################################

immutable CharacterIterator{T <: AbstractString}
    s::T
end
Base.start(x::CharacterIterator) = start(x.s)
Base.next(x::CharacterIterator, i::Integer) = next(x.s, i)
Base.done(x::CharacterIterator, i::Integer) = done(x.s, i)
Base.length(x::CharacterIterator) = length(x.s)

Base.nextind(x::CharacterIterator, i::Integer) = nextind(x.s, i)
Base.chr2ind(x::CharacterIterator, i::Integer) = chr2ind(x.s, i::Integer)
iteratortype{T <: CharacterIterator}(::Type{T}) = CharacterIterator
Base.convert{T}(::Type{CharacterIterator{T}}, x::T) = CharacterIterator(x)

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
iteratortype{T <: GraphemeIterator}(::Type{T}) = GraphemeIterator
Base.convert{T}(::Type{GraphemeIterator{T}}, x::T) = GraphemeIterator(x)


typealias StringIterator{T} Union{GraphemeIterator{T}, CharacterIterator{T}}
Base.endof(g::StringIterator) = endof(g.s)
Base.SubString(x::StringIterator, i, j) = iteratortype(x)(SubString(x.s, i, j))
Base.string(x::StringIterator) = x.s
Base.isempty(x::StringIterator) = isempty(x.s)
Base.eltype{T}(x::StringIterator{T}) = T
Base.isless(x::StringIterator, y::StringIterator) = isless(x.s, y.s)
Base.search(x::StringIterator, args...) = search(x.s, args...)
Base.searchsortedfirst(x::StringIterator, args...) = searchsortedfirst(x.s, args...)
Base.searchsortedlast(x::StringIterator, args...) = searchsortedlast(x.s, args...)
Base.searchsorted(x::StringIterator, args...) = searchsorted(x.s, args...)
iteratortype(x::StringIterator) = iteratortype(typeof(x))

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
        function $x(dist::PreMetric, s1::AbstractString, s2::AbstractString)
            $x(dist, CharacterIterator(s1), CharacterIterator(s2))
        end

        function $x(dist::PreMetric, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
            $x(dist, CharacterIterator(s1), CharacterIterator(s2), len1, len2)
        end

        function $x(dist::PreMetric, s1::StringIterator, s2::StringIterator)
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

function compare(dist::PreMetric, s1::StringIterator, s2::StringIterator, 
    len1::Integer, len2::Integer)
    1.0 - evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::StringIterator, s2::StringIterator,
    len1::Integer, len2::Integer)
    distance = evaluate(dist, s1, s2, len1, len2)
    len2 == 0 ? 1.0 : 1.0 - distance / len2
end

# compare always return a value between 0 and 1. 
# When string length < q for qgram distance, returns s1 == s2
function compare(dist::AbstractQGram, 
    s1::StringIterator, s2::StringIterator, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    evaluate(dist, s1, s2, len1, len2)
end

function compare(dist::QGram, 
    s1::StringIterator, s2::StringIterator, 
    len1::Integer, len2::Integer)
    len1 <= (dist.q - 1) && return convert(Float64, s1 == s2)
    distance = evaluate(dist, s1, s2, len1, len2)
    1 - distance / (len1 + len2 - 2 * dist.q + 2)
end





end 