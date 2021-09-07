
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
# Not really needed but avoid multi-encapsulation
string_with_length(s::StringWithLength) = s
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
    l = 0
    for (ch1, ch2) in zip(s1, s2)
        ch1 != ch2 && break
        l += 1
    end
    return l
end


