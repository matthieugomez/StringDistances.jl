# Return start of commn substring in s1, start of common substring in s2, and length of substring
# Indexes refer to character number, not index (differ for Unicode strings)
function longest_common_substring(s1::AbstractString, s2::AbstractString)
    if length(s1) > length(s2)
        start2, start1, len = longest_common_substring(s2, s1)
    else
        start1, start2, len = 0, 0, 0
        p = zeros(Int, length(s2))
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

function matching_blocks!(x::Set{Tuple{Int, Int, Int}}, s1::AbstractString, s2::AbstractString, start1::Integer, start2::Integer)
    a = longest_common_substring(s1, s2)
    if a[3] > 0
        push!(x, (a[1] + start1 - 1, a[2] + start2 - 1, a[3]))
        s1before = SubString(s1, firstindex(s1), nextind(s1, 0, a[1] - 1))
        s2before = SubString(s2, firstindex(s2), nextind(s2, 0, a[2] - 1))
        matching_blocks!(x, s1before, s2before, start1, start2)
        if ((a[1] + a[3]) <= lastindex(s1)) & ((a[2] + a[3]) <= lastindex(s2))
            s1after = SubString(s1, nextind(s1, 0, a[1] + a[3]), lastindex(s1))
            s2after = SubString(s2, nextind(s2, 0, a[2] + a[3]), lastindex(s2))
            matching_blocks!(x, s1after, s2after, start1 + a[1] + a[3] - 1, start2 + a[2] + a[3] - 1)
        end
    end
end

function matching_blocks(s1::AbstractString, s2::AbstractString)
    x = Set{Tuple{Int, Int, Int}}()
    matching_blocks!(x, s1, s2, 1, 1)
    return x
end

struct RatcliffObershelp <: PreMetric end

function evaluate(dist::RatcliffObershelp, s1::AbstractString, s2::AbstractString)
    n_matched = 0
    for x in matching_blocks(s1, s2)
        n_matched += x[3]
    end
    1.0 - 2 * n_matched / (length(s1) + length(s2))
end
