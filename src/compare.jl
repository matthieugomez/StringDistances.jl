##############################################################################
##
## compare
## compare always return a value between 0 and 1.
##
##############################################################################

function compare(dist::PreMetric, s1::AbstractString, s2::AbstractString)
    1.0 - evaluate(dist, s1, s2)
end

function compare(dist::Union{Hamming, Levenshtein, DamerauLevenshtein}, 
    s1::AbstractString, s2::AbstractString)
    len = max(length(s1), length(s2))
    len == 0 ? 1.0 : 1.0 - evaluate(dist, s1, s2) / len
end

function compare(dist::AbstractQGram, s1::AbstractString, s2::AbstractString)
    # When string length < q for qgram distance, returns s1 == s2
    len1 = length(s1) ; len2 = length(s2)
    min(len1, len2) <= (dist.q - 1) && return convert(Float64, s1 == s2)
    if typeof(dist) <: QGram
        1 - evaluate(dist, s1, s2) / (len1 + len2 - 2 * dist.q + 2)
    else
        1 - evaluate(dist, s1, s2)
    end
end
