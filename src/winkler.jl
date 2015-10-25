##############################################################################
##
## Winkler
##
##############################################################################

type Winkler{T1, T2 <: Real, T3 <: Real}
    dist::T1
    scaling_factor::T2      # scaling factor. Default to 0.1
    boosting_limit::T3      # boost threshold. Default to 1.0 
end

# restrict to distance between 0 and 1
Winkler(x) = Winkler(x, 0.1, 1.0)

function evaluate(dist::Winkler, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    distance = evaluate(Normalized(dist.dist), s1, s2, len1, len2)
    l = common_prefix(s1, s2, 4)[1]
    # common prefix adjustment
    if distance <= dist.boosting_limit
        distance -= distance * l * dist.scaling_factor
    end
    return distance
end