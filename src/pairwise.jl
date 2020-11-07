_allocmatrix(X, Y, T) = Matrix{T}(undef, length(X), length(Y))
_allocmatrix(X, T) = Matrix{T}(undef, length(X), length(X))

import Distances: pairwise

@doc """
    pairwise(dist::StringDistance, itr; eltype = Float64, precalc = nothing)
    pairwise(dist::StringDistance, itr1, itr2; eltype = Float64, precalc = nothing)

`pairwise` returns the distance matrix between all pairs of elements in `itr`
according to the distance `dist`. The element type of the returned matrix
can be set via `eltype`. For QGramDistances precalculation will be used either
if `precalc` is set to true or if there are more than 5 elements in `itr`.
Set `precalc` to false if no precalculation should be used, regardless of length.

Both symmetric and asymmetric versions are available.

### Examples
```julia-repl
julia> using StringDistances
julia> iter = ["New York", "Princeton"]
julia> pairwise(Levenshtein(), iter) # symmetric
2×2 Array{Float64,2}:
 0.0  9.0
 9.0  0.0
julia> iter2 = ["San Francisco"]
julia> pairwise(Levenshtein(), iter, iter2) # asymmetric
2×1 Array{Float64,2}:
 12.0
 10.0
```
"""
pairwise

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
