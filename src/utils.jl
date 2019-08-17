function reorder(s1::AbstractString, s2::AbstractString)
    len1 = length(s1)
    len2 = length(s2)
    if len2 > len1
        return s2, len2, s1, len1
    else
        return s1, len1, s2, len2
    end
end

## Find common prefixes (up to lim. -1 means Inf)
function common_prefix(s1::AbstractString, s2::AbstractString, lim::Integer = -1)
    x1 = iterate(s1)
    x2 = iterate(s2)
    l = 0
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
