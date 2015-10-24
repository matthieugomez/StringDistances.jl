##############################################################################
##
## Gram Iterator iterates through q-grams of a string
##
##############################################################################

type QGramIterator{S <: AbstractString, T <: Integer}
    s::S
    q::T
end
function Base.start(qgram::QGramIterator)
    len = length(qgram.s)
    (1, len < qgram.q ? endof(qgram.s) + 1 : chr2ind(qgram.s, qgram.q))
end
function Base.next(qgram::QGramIterator, state)
    istart, iend = state
    element = convert(typeof(qgram.s), SubString(qgram.s, istart, iend))
    nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
    return element, nextstate
end
function Base.done(qgram::QGramIterator, state)
    istart, idend = state
    done(qgram.s, idend)
end
Base.eltype(qgram::QGramIterator) = typeof(qgram.s)
Base.length(qgram::QGramIterator) = length(qgram.s) - qgram.q + 1

##############################################################################
##
## A Bag is a Set that allows duplicated values
## Implemented as Dictionary from elements => number of duplicates
##
##############################################################################

type Bag{Tv, Ti <: Integer}
    dict::Dict{Tv, Ti}
    Bag() = new(Dict{Tv, Ti}())
end

function Base.push!{Tv, Ti}(bag::Bag{Tv, Ti}, x::Tv)
    bag.dict[x] = get(bag.dict, x, zero(Ti)) + one(Ti)
    return bag
end

function Base.delete!{Tv, Ti}(bag::Bag{Tv, Ti}, x::Tv)
    v = get(bag.dict, x, zero(Ti))
    if v > zero(Ti)
        bag.dict[x] = v - one(Ti)
    end
    return x
end

Base.length(bag::Bag) = convert(Int, sum(values(bag.dict)))

function Bag(s)
    bag = Bag{eltype(s), UInt}()
    for x in s
        push!(bag, x)
    end
    return bag
end

##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

type QGram{T <: Integer}
    q::T
end
QGram() = QGram(2)

function evaluate(dist::QGram, s1::AbstractString, s2::AbstractString) 
    len1, len2 = length(s1), length(s2)
    len1 > len2 && return evaluate(dist, s2, s1)
    len2 == 0 && return 0

    q1 = QGramIterator(s1, dist.q)
    q2 = QGramIterator(s2, dist.q)
    bag = Bag(q2)
    for ch in q1
        delete!(bag, ch)
    end
    # number non matched in s1 : n1 - (n2 - length(bag)) 
    # number non matched in s2 : length(bag)
    return length(q1) - length(q2) + 2 * length(bag)
end

qgram(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(QGram(q), s1::AbstractString, s2::AbstractString)

##############################################################################
##
## cosine 
##
## 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
##############################################################################

type Cosine{T <: Integer}
    q::T
end
Cosine() = Cosine(2)

function evaluate(dist::Cosine, s1::AbstractString, s2::AbstractString) 
    len1, len2 = length(s1), length(s2)
    len1 > len2 && return evaluate(dist, s2, s1)
    len2 == 0 && return 0.0

    bag2 = Bag(QGramIterator(s2, dist.q))
    bag1 = Bag(QGramIterator(s1, dist.q))
    numerator = 0
    for (k, v1) in bag1.dict
        numerator += v1 * get(bag2.dict, k, 0)
    end
    denominator = sqrt(sumabs2(values(bag1.dict))) * sqrt(sumabs2(values(bag2.dict)))
    denominator == 0 ? 1.0 : 1.0 - numerator / denominator
end

cosine(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(Cosine(q), s1::AbstractString, s2::AbstractString)

##############################################################################
##
## Jaccard
##
## Denote Q(s, q) the set of tuple of length q in s
## 1 - |intersect(Q(s1, q), Q(s2, q))| / |union(Q(s1, q), Q(s2, q))|
##
##############################################################################

type Jaccard{T <: Integer}
    q::T
end
Jaccard() = Jaccard(2)

function evaluate(dist::Jaccard, s1::AbstractString, s2::AbstractString) 
    len1, len2 = length(s1), length(s2)
    len1 > len2 && return evaluate(dist, s2, s1)
    len2 == 0 && return 0.0

    set2 = Set(QGramIterator(s2, dist.q))
    set1 = Set(QGramIterator(s1, dist.q))
    numerator = 0
    for x in set1
        if x in set2
            numerator += 1
        end
    end
    denominator = length(set1) + length(set2) - numerator
    return 1.0 - numerator / denominator
end

jaccard(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(Jaccard(q), s1::AbstractString, s2::AbstractString)