##############################################################################
##
## TokenSort
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
type TokenSort{T <: PreMetric} <: PreMetric
    dist::T
end

function compare{T <: AbstractString}(dist::TokenSort, s1::T, s2::T, len1::Integer, len2::Integer)
    s1 = join(sort!(split(s1)), " ")
    s2 = join(sort!(split(s2)), " ")
    compare(dist.dist, s1, s2)
end

##############################################################################
##
## TokenSet
## http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/
##
##############################################################################
type TokenSet{T <: PreMetric} <: PreMetric
    dist::T
end

function compare{T <: AbstractString}(dist::TokenSet, s1::T, s2::T, len1::Integer, len2::Integer)
    v0, v1, v2 = _separate!(split(s1), split(s2))
    s0 = join(v0, " ")
    s1 = join(chain(v0, v1), " ")
    s2 = join(chain(v0, v2), " ")
    if isempty(s0)
        # otherwise compare(dist, "", "a")== 1.0 
        compare(dist.dist, s1, s2)
    else
        max(compare(dist.dist, s0, s1), 
            compare(dist.dist, s1, s2), 
            compare(dist.dist, s0, s2))        
    end
end

# separate 2 vectors in intersection, setdiff1, setdiff2 (all sorted)
function _separate!(v1::Vector, v2::Vector)
    sort!(v1)
    sort!(v2)
    out = eltype(v1)[]
    start = 1
    i1 = 0
    while i1 < length(v1)
        i1 += 1
        x = v1[i1]
        i2 = searchsortedfirst(v2, x, start, length(v2), Base.Forward)
        i2 > length(v2) && break 
        if i2 > 0 && v2[i2] == x
            deleteat!(v1, i1)
            deleteat!(v2, i2)
            push!(out, x)
            i1 -= 1
            start = i2 
        end
    end
    return out, v1, v2
end
