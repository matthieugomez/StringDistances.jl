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
TokenMax,
graphemes2

##############################################################################
##
## Iterator
##
##############################################################################

isgraphemebreak(c1::Char, c2::Char) =
    ccall(:utf8proc_grapheme_break, Bool, (UInt32, UInt32), c1, c2)

immutable GraphemeIterator2{S<:AbstractString}
    s::S # original string (for generation of SubStrings)
end
graphemes2(s::AbstractString) = GraphemeIterator2{typeof(s)}(s)

Base.eltype{S}(::Type{GraphemeIterator2{S}}) = SubString{S}

function Base.length(g::GraphemeIterator2)
    c0 = Char(0x00ad) # soft hyphen (grapheme break always allowed after this)
    n = 0
    for c in g.s
        n += isgraphemebreak(c0, c)
        c0 = c
    end
    return n
end

Base.start(g::GraphemeIterator2) = start(g.s)
Base.done(g::GraphemeIterator2, i) = done(g.s, i)

function Base.next(g::GraphemeIterator2, i)
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

# functions not defined in base
Base.nextind(g::GraphemeIterator2, state::Integer) = next(g, state)[2]
function Base.chr2ind(g::GraphemeIterator2, idx::Integer)
    state = start(g)
    i = 0
    while !done(g, state)
        i += 1
        i == idx && return state
        ch, state = next(g, state)
    end
end
Base.endof(g::GraphemeIterator2) = endof(g.s)

typealias GraphemeOrString Union{GraphemeIterator2, AbstractString}
Base.SubString(x::GraphemeIterator2, i, j) = SubString(x.s, i, j)
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