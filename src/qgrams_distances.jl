
##############################################################################
##
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## cosine is 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##############################################################################


##############################################################################
##
## A Bag is like Set that it allows duplicated values
## I implement it as dictionary from elements => number of duplicates
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

function Bag(s::AbstractString, q::Integer)
    bag = Bag{typeof(s), UInt}()
    @inbounds for i in 1:(length(s) - q + 1)
        push!(bag, s[i:(i + q - 1)])
    end
    return bag
end

##############################################################################
##
## q-gram 
##
##############################################################################

type QGram{T <: Integer}
    q::T
end

function evaluate{T}(dist::QGram, s1::T, s2::T) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0
    n2 = length(s2) - dist.q + 1
    bag = Bag(s2, dist.q)
    count = 0
    n1 = length(s1) - dist.q + 1
    for i1 in 1:n1
        @inbounds ch = s1[i1:(i1 + dist.q - 1)]
        delete!(bag, ch)
    end
    # number non matched in s1 : n1 - (n2 - length(bag)) 
    # number non matched in s2 : length(bag)
    return n1 - n2 + 2 * length(bag)
end

qgram{T}(s1::T, s2::T; q = 2) = evaluate(QGram(q), s1, s2)

##############################################################################
##
## cosine 
##
##############################################################################

type Cosine{T <: Integer}
    q::T
end

function evaluate{T}(dist::Cosine, s1::T, s2::T) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0.0

    bag2 = Bag(s2, dist.q)
    bag1 = Bag(s1, dist.q)

    count = 0
    for (k, v1) in bag1.dict
        count += v1 * get(bag2.dict, k, 0)
    end
    denominator = sqrt(sumabs2(values(bag1.dict))) * sqrt(sumabs2(values(bag2.dict)))
    denominator == 0 ? 1.0 : 1.0 - count / denominator
end

cosine{T}(s1::T, s2::T; q = 2) = evaluate(Cosine(q), s1, s2)

##############################################################################
##
## Jaccard
##
## Denote Q(s, q) the set of tuple of length q in s
## jaccard(s1, s2, q) = 1 - |intersect(Q(s1, q), Q(s2, q))| / |union(Q(s1, q), Q(s2, q))|
##
##############################################################################

type Jaccard{T <: Integer}
    q::T
end

function evaluate{T}(dist::Jaccard, s1::T, s2::T) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0.0

   
    set2 = Set{T}()
    n2 = length(s2) - dist.q + 1
    @inbounds for i2 in 1:n2
        push!(set2, s2[i2:(i2 + dist.q - 1)])
    end

    set1 = Set{T}()
    n1 = length(s1) - dist.q + 1
    @inbounds for i1 in 1:n1
        push!(set1, s1[i1:(i1 + dist.q - 1)])
    end

    n_intersect = 0
    for x in set1
        if x in set2
            n_intersect += 1
        end
    end

    return 1.0 - n_intersect / (length(set1) + length(set2) - n_intersect)
end

jaccard{T}(s1::T, s2::T; q = 2) = evaluate(Jaccard(q), s1, s2)




