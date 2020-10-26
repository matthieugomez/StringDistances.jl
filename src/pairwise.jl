_allocmat(X, Y, T) = Matrix{T}(undef, length(X), length(Y))
_allocmat(X, T) = Matrix{T}(undef, length(X), length(X))

pairwise(dist::StringDistance, X, Y; eltype = Float64) = 
    pairwise!(_allocmat(X, Y, eltype), dist, X, Y)

pairwise(dist::StringDistance, X; eltype = Float64) = 
    pairwise!(_allocmat(X, eltype), dist, X)

function pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X) where {N<:Number}
    if dist isa SemiMetric
        _symmetric_pairwise!(R, dist, X)
    else
        _asymmetric_pairwise!(R, dist, X, X)
    end
end

function pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X, Y) where {N<:Number}
    _asymmetric_pairwise!(R, dist, X, Y)
end

_precalc(X, PT, q) = PT[PT(X[i], q) for i in 1:length(X)]

const PrecalcMinLength = 5 # Only precalc if length >= 5

function _symmetric_pairwise!(R, dist::QGramDistance, X; precalc = nothing, precalcType = QGramSortedVector)
    # precalc if set to true or if isnothing and length is at least min length
    shouldprecalc = (precalc === true) | (isnothing(precalc) & length(X) >= PrecalcMinLength)
    objs = shouldprecalc ? _precalc(X, precalcType, q(dist)) : X

    for i in 1:length(objs)
        R[i, i] = 0
        for j in (i+1):length(objs)
            R[i, j] = R[j, i] = evaluate(dist, objs[i], objs[j])
        end
    end
    return R
end
