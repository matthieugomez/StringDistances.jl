
##############################################################################
##
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## cosine is 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##############################################################################


##############################################################################
##
## A bag is a Set with repeated values
##
##############################################################################

type Bag{Tv <: Union{Char, AbstractString}, Ti <: Integer}
    dict::Dict{Tv, Ti}
end

function Bag(s::AbstractString, q::Integer)
    dict = Dict{typeof(s), Int}()
    for i in 1:(length(s) - q + 1)
        ch = s[i:(i + q - 1)]
        dict[ch] = get(dict, ch, 0) + 1
    end
    return Bag(dict)
end

Base.in{Tv <: Union{Char, AbstractString}, Ti}(x::Tv, bag::Bag{Tv, Ti}) = get(bag.dict, x, 0) > 0
function Base.pop!{Tv, Ti}(bag::Bag{Tv, Ti}, x::Tv)
    bag.dict[x] -= 1
    return x
end
Base.length(bag::Bag) = sum(values(bag.dict))

##############################################################################
##
## q-gram 
##
##############################################################################

type QGram{T <: Integer}
    q::T
end

function evaluate(dist::QGram, s1::AbstractString, s2::AbstractString) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0
    
    # set with repeated element : for each key, number of repetitions
    bag = Bag(s2, dist.q)
    count = 0
    n1 = length(s1) - dist.q + 1
    for i1 in 1:n1
        ch = s1[i1:(i1 + dist.q - 1)]
        if ch in bag
            pop!(bag, ch)
            count += 1
        end
    end

    return n1 - count + length(bag)
end

qgram(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(QGram(q), s1, s2)

##############################################################################
##
## cosine 
##
##############################################################################

type Cosine{T <: Integer}
    q::T
end

function evaluate(dist::Cosine, s1::AbstractString, s2::AbstractString) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0.0

    bag2 = Bag(s2, dist.q)
    bag1 = Bag(s1, dist.q)

    count = 0
    for x in keys(bag1.dict)
        if x in bag2
            count += bag1.dict[x] * bag2.dict[x]
        end
    end

    return 1.0 - count / (sqrt(sumabs2(values(bag1.dict))) * sqrt(sumabs2(values(bag2.dict))))
end

cosine(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(Cosine(q), s1, s2)

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

function evaluate(dist::Jaccard, s1::AbstractString, s2::AbstractString) 
    length(s1) > length(s2) && return evaluate(dist, s2, s1)
    length(s2) == 0 && return 0.0

    n2 = length(s2) - dist.q + 1
    n1 = length(s1) - dist.q + 1
   
    set2 = Set{typeof(s2)}()
    for i2 in 1:n2
        push!(set2, s2[i2:(i2 + dist.q - 1)])
    end

    set1 = Set{typeof(s1)}()
    for i1 in 1:n1
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

jaccard(s1::AbstractString, s2::AbstractString; q = 2) = evaluate(Jaccard(q), s1, s2)


