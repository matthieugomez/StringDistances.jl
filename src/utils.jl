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


# Return start of commn substring in s1, start of common substring in s2, and length of substring
# Indexes refer to character number, not index
function longest_common_substring(s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    if len1 > len2
        start2, start1, len = longest_common_substring(s2, s1, len2, len1)
    else
        start1, start2, len = 0, 0, 0
        p = zeros(Int, len2)
        i1 = 0
        for ch1 in s1
            i1 += 1
            i2 = 0
            oldp = 0
            for ch2 in s2
                i2 += 1
                newp = 0
                if ch1 == ch2
                    newp = oldp > 0 ? oldp : i2
                    currentlength = (i2 - newp + 1)
                    if currentlength > len
                        start1, start2, len = i1 - currentlength + 1, newp, currentlength
                    end
                end
                p[i2], oldp = newp, p[i2]
            end
        end
    end
    return start1, start2, len
end
