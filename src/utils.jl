# String with Length
# This allows to compute length once and only once
struct StringWithLength{T<:AbstractString} <: AbstractString
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
    if length(s1) <= length(s2)
         return s1, s2
    else
        return s2, s1
    end
end

 
## Find common prefixes (up to lim. -1 means Inf)
function remove_prefix(s1::AbstractString, s2::AbstractString, lim::Integer = -1)
    l = 0
    x1 = iterate(s1)
    x2 = iterate(s2)
    while (x1 !== nothing) & (x2 !== nothing) & (l < lim || lim < 0)
        ch1, state1 = x1
        ch2, state2 = x2
        ch1 != ch2 && break
        x1 = iterate(s1, state1)
        x2 = iterate(s2, state2)
        l += 1
    end
    return l, x1, x2
end