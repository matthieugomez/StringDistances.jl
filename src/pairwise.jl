_allocmatrix(X, Y, T) = Matrix{T}(undef, length(X), length(Y))
_allocmatrix(X, T) = Matrix{T}(undef, length(X), length(X))

@doc """
    pairwise(dist::StringDistance, itr; eltype = Float64, preprocess = nothing)
    pairwise(dist::StringDistance, itr1, itr2; eltype = Float64, preprocess = nothing)

Compute distances between all pairs of elements in `itr`according to the `StringDistance` 
`dist`. The element type of the returned distance matrix can be set via `eltype`. 

For QGramDistances preprocessing will be used either if `preprocess` is set to true or 
if there are more than 5 elements in `itr`. Set `preprocess` to false if no 
preprocessing should be used, regardless of length.

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
Distances.pairwise

Distances.pairwise(dist::StringDistance, X, Y; eltype = Float64, preprocess = nothing) =
    pairwise!(_allocmatrix(X, Y, eltype), dist, X, Y; preprocess = preprocess)

Distances.pairwise(dist::StringDistance, X; eltype = Float64, preprocess = nothing) =
    pairwise!(_allocmatrix(X, eltype), dist, X; preprocess = preprocess)

pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X; preprocess = nothing) where {N<:Number} =
    (dist isa SemiMetric) ?
        _symmetric_pairwise!(R, dist, X; preprocess = preprocess) :
        _asymmetric_pairwise!(R, dist, X, X; preprocess = preprocess)

pairwise!(R::AbstractMatrix{N}, dist::StringDistance, X, Y; preprocess = nothing) where {N<:Number} =
    _asymmetric_pairwise!(R, dist, X, Y; preprocess = preprocess)

_preprocess(X, PT, q) = PT[PT(X[i], q) for i in 1:length(X)]

const PrecalcMinLength = 5 # Only precalc if length >= 5

preprocess_if_needed(X, dist::StringDistance, preprocess, preprocessType) = X

function preprocess_if_needed(X, dist::QGramDistance, preprocess, preprocessType)
    # preprocess only if a QGramDistance and
    # if precalc set to true or if isnothing and length is at least min length
    cond = (preprocess === true) ||
                (isnothing(preprocess) && length(X) >= PrecalcMinLength)
    cond ? _preprocess(X, preprocessType, dist.q) : X
end

function _symmetric_pairwise!(R, dist::StringDistance, X;
    preprocess = nothing, preprocessType = QGramSortedVector)

    objs = preprocess_if_needed(X, dist, preprocess, preprocessType)

    for i in 1:length(objs)
        R[i, i] = 0
        Threads.@threads for j in (i+1):length(objs)
            R[i, j] = R[j, i] = evaluate(dist, objs[i], objs[j])
        end
    end
    return R
end

function _asymmetric_pairwise!(R, dist::StringDistance, X, Y;
    preprocess = nothing, preprocessType = QGramSortedVector)

    objsX = preprocess_if_needed(X, dist, preprocess, preprocessType)
    objsY = preprocess_if_needed(Y, dist, preprocess, preprocessType)

    for i in 1:length(objsX)
        Threads.@threads for j in 1:length(objsY)
            R[i, j] = evaluate(dist, objsX[i], objsY[j])
        end
    end
    return R
end
