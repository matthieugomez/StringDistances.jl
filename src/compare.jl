##############################################################################
##
## compare
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    1.0 - evaluate(dist, s1, s2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::AbstractString, s2::AbstractString)
    distance = evaluate(dist, s1, s2)
    len = max(length(s1), length(s2))
    len == 0 ? 1.0 : 1.0 - distance / len
end

# compare always return a value between 0 and 1.
# When string length < q for qgram distance, returns s1 == s2
function compare(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
    len = min(length(s1), length(s2))
    len <= (dist.q - 1) && return convert(Float64, s1 == s2)
    evaluate(dist, s1, s2)
end

function compare(dist::QGram, s1::AbstractString, s2::AbstractString)
    len1 = length(s1)
    len2 = length(s2)
    min(len1, len2) <= (dist.q - 1) && return convert(Float64, s1 == s2)
    distance = evaluate(dist, s1, s2)
    1 - distance / (len1 + len2 - 2 * dist.q + 2)
end
