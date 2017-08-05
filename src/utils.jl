function reorder(s1::AbstractString, s2::AbstractString)
    len1 = length(s1)
    len2 = length(s2)
    if len2 > len1
        return s2, len2, s1, len1
    else
        return s1, len1, s2, len2
    end
end
