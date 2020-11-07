_allocmatrix(X, Y, T) = Matrix{T}(undef, length(X), length(Y))
_allocmatrix(X, T) = Matrix{T}(undef, length(X), length(X))

pairwise(dist::StringDistance, X, Y; eltype = Float64, precalc = nothing) =
    pairwise!(_allocmatrix(X, Y, eltype), dist, X, Y; precalc)

pairwise(dist::StringDistance, X; eltype = Float64, precalc = nothing) =
    pairwise!(_allocmatrix(X, eltype), dist, X; precalc)

pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X; precalc = nothing) where {N<:Number} =
    (dist isa SemiMetric) ?
        _symmetric_pairwise!(R, dist, X; precalc) :
        _asymmetric_pairwise!(R, dist, X, X; precalc)

pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X, Y; precalc = nothing) where {N<:Number} =
    _asymmetric_pairwise!(R, dist, X, Y; precalc)

_precalc(X, PT, q) = PT[PT(X[i], q) for i in 1:length(X)]

const PrecalcMinLength = 5 # Only precalc if length >= 5

function precalc_if_needed(X, dist::StringDistance, precalc, precalcType)
    # precalc only if a QGramDistance and
    # if precalc set to true or if isnothing and length is at least min length
    !isa(dist, QGramDistance) && return X
    cond = (precalc === true) ||
                (isnothing(precalc) & length(X) >= PrecalcMinLength)
    cond ? _precalc(X, precalcType, dist.q) : X
end

function _symmetric_pairwise!(R, dist::StringDistance, X;
    precalc = nothing, precalcType = QGramSortedVector)

    objs = precalc_if_needed(X, dist, precalc, precalcType)

    for i in 1:length(objs)
        R[i, i] = 0
        Threads.@threads for j in (i+1):length(objs)
            R[i, j] = R[j, i] = evaluate(dist, objs[i], objs[j])
        end
    end
    return R
end

function _asymmetric_pairwise!(R, dist::StringDistance, X, Y;
    precalc = nothing, precalcType = QGramSortedVector)

    objsX = precalc_if_needed(X, dist, precalc, precalcType)
    objsY = precalc_if_needed(Y, dist, precalc, precalcType)

    for i in 1:length(objsX)
        Threads.@threads for j in 1:length(objsY)
            R[i, j] = evaluate(dist, objsX[i], objsY[j])
        end
    end
    return R
end
