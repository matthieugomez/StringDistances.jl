##############################################################################
##
## Winkler
##
##############################################################################

struct Winkler{T1 <: PreMetric, T2 <: Real, T3 <: Real} <: PreMetric
    dist::T1
    scaling_factor::T2      # scaling factor. Default to 0.1
    boosting_limit::T3      # boost threshold. Default to 0.7
end

# restrict to distance between 0 and 1
Winkler(x) = Winkler(x, 0.1, 0.7)

function compare(dist::Winkler, s1::AbstractStringorGraphemeIterator, s2::AbstractStringorGraphemeIterator, len1::Integer, len2::Integer)
    score = compare(dist.dist, s1, s2, len1, len2)
    l = common_prefix(s1, s2, 4)[1]
    # common prefix adjustment
    if score >= dist.boosting_limit
        score += l * dist.scaling_factor * (1 - score)
    end
    return score
end