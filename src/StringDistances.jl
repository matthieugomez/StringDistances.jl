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
TokenMax,
graphemeiterator
##############################################################################
##
## Define GraphemeIterator as AbstractString
##
## Argument for non inheritance: AbstractString is generally something that gives char as individual. 
## Argument for inheritance: prevind, nextind, chr2ind are defined once start, next, done, isvalid, endof are defined. Also this allows to define functions with AbstractString signature & having SubString work.
##############################################################################
# from Base. I redefine it because I want AbstractStringinheritance
immutable GraphemeIterator{S<:AbstractString} <: AbstractString
    s::S # original string (for generation of SubStrings)
end
graphemeiterator(s::AbstractString) = GraphemeIterator{typeof(s)}(s)
eltype{S}(::Type{GraphemeIterator{S}}) = SubString{S}
function length(g::GraphemeIterator)
    c0 = Char(0x00ad) # soft hyphen (grapheme break always allowed after this)
    n = 0
    for c in g.s
        n += isgraphemebreak(c0, c)
        c0 = c
    end
    return n
end
start(g::GraphemeIterator) = start(g.s)
done(g::GraphemeIterator, i::Int) = done(g.s, i)
function next(g::GraphemeIterator, i::Int)
    s = g.s
    j = i
    c0, k = next(s, i)
    while !done(s, k) # loop until next grapheme is s[i:j]
        c, ℓ = next(s, k)
        isgraphemebreak(c0, c) && break
        j = k
        k = ℓ
        c0 = c
    end
    return (SubString(s, i, j), k)
end
==(g1::GraphemeIterator, g2::GraphemeIterator) = g1.s == g2.s
hash(g::GraphemeIterator, h::UInt) = hash(g.s, h)
isless(g1::GraphemeIterator, g2::GraphemeIterator) = isless(g1.s, g2.s)
show{S}(io::IO, g::GraphemeIterator{S}) = print(io, "length-$(length(g)) GraphemeIterator{$S} for \"$(g.s)\"")


# added
#used in prevind nextind
function Base.isvalid(s::GraphemeIterator, i::Integer)
    if !isvalid(s.s, i) 
        return false
    else
        i0 = prevind(s.s, i)
        return i0 < start(s.s) || isgraphemebreak(s.s[i0], s.s[i])
    end
end
function endof(s::GraphemeIterator)
    c0 = Char(0x00ad)
    i = endof(s.s)
    i0 = start(s.s)
    while i >= i0 && !isgraphemebreak(s.s[i], c0)
        i = prevind(s.s, i)
        c0 = s.s[i]
    end
    i
end

# 1. issues with stuff like search  or print_escaped where character iteration vs string iteration matters. I need to pass the original string for now
Base.search(x::GraphemeIterator, s::Vector{Char}) = search(x.s, s)
# 2. issue with keeping iterator property for stuff like split, join. for now, I decide to loose the enumerator property but add it back after join. But SubString for instance does not loose the property
Base.split(x::GraphemeIterator, args...) = split(x.s, args...)
iterator{T <: GraphemeIterator}(::Type{T}, x::AbstractString) = graphemeiterator(x)
iterator{T <: AbstractString}(::Type{T}, x::AbstractString) = x

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