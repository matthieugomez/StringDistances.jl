function evaluate(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    len1, len2 = length(s1), length(s2)
    if len1 > len2
        return evaluate(dist, s2, s1, len2, len1)
    else
        return evaluate(dist, s1, s2, len1, len2)
    end
end
