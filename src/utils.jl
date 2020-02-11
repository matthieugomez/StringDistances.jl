
##############################################################################
##
## Some things about Strings

# length: number of characters
# ncodeunits: Return the number of code units in a string (aking to index of vector). 
# Not all such indices  are valid – they may not be the start of a character,.
# sizeof:  Size, in bytes, of the string str. Equal to the number of code units in str  
# multiplied by the size, in bytes, of one code unit in str.

# lastindex: Return the last index of a collection
# nextinds(s, i):  return the index of the start of the character whose encoding starts after index i
# nextind(s, 0, N): return the index of the Nth character of s (or, if there are 
# less than N characters, return ncodeunits(str) + (N - length(s))

##############################################################################


# This type allows to compute length once and for all
struct StringWithLength{T <: AbstractString} <: AbstractString
    s::T
    l::Int
end
string_with_length(s::AbstractString) = StringWithLength(s, length(s))
Base.length(s::StringWithLength) = s.l
Base.iterate(s::StringWithLength, i::Integer = firstindex(s.s)) = iterate(s.s, i)
Base.nextind(s::StringWithLength, i::Int, n::Int = 1) = nextind(s.s, i, n)
Base.ncodeunits(s::StringWithLength) = ncodeunits(s.s)
Base.isvalid(s::StringWithLength, i::Int) = isvalid(s.s, i)


function reorder(s1::AbstractString, s2::AbstractString)
    s1 = string_with_length(s1)
    s2 = string_with_length(s2)
    (length(s1) <= length(s2)) ? (s1, s2) : (s2, s1)
end

function reorder(s1, s2)
    (length(s1) <= length(s2)) ? (s1, s2) : (s2, s1)
end

function common_prefix(s1, s2)
    x1 = iterate(s1)
    x2 = iterate(s2)
    l = 0
    while (x1 !== nothing) & (x2 !== nothing)
        ch1, state1 = x1
        ch2, state2 = x2
        ch1 != ch2 && break
        l += 1
        x1 = iterate(s1, state1)
        x2 = iterate(s2, state2)
    end
    return l, x1, x2
end



function _take(s, n::Integer)
    Base.Iterators.take(s, n)
end
function _take(s::AbstractString, n::Integer)
   SubString(s, firstindex(s), nextind(s, 0, n))
end

function _drop(s, n::Integer)
    Base.Iterators.drop(s, n)
end
function _drop(s::AbstractString, n::Integer)
   SubString(s, nextind(s, 0, n + 1), lastindex(s))
end

function _slice(s, n1::Integer, n2::Integer)
    Base.Iterators.take(Base.Iterators.drop(s, n1), n2 - n1)
end
function _slice(s::AbstractString, n1::Integer, n2::Integer)
   SubString(s, nextind(s, 0, n1 + 1),  nextind(s, 0, n2))
end


