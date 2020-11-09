@doc """
    pairwise(dist::StringDistance, itr; preprocess = nothing)
    pairwise(dist::StringDistance, itr1, itr2; preprocess = nothing)

Compute distances between all pairs of elements in `itr` according to the
`StringDistance` `dist`.

For QGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `itr`. Set `preprocess` to 
false if no preprocessing should be used, regardless of length.

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

function Distances.pairwise(dist::StringDistance, X, Y; preprocess = nothing)
    T = result_type(dist, eltype(X), eltype(Y))
    R = Matrix{T}(undef, length(X), length(Y))
    pairwise!(R, dist, X, Y; preprocess = preprocess)
end

function Distances.pairwise(dist::StringDistance, X; preprocess = nothing)
    T = result_type(dist, eltype(X), eltype(X))
    R = Matrix{T}(undef, length(X), length(X))
    pairwise!(R, dist, X; preprocess = preprocess)
end

@doc """
    pairwise!(r::AbstractMatrix, dist::StringDistance, itr; preprocess = nothing)
    pairwise!(r::AbstractMatrix, dist::StringDistance, itr1, itr2; preprocess = nothing)

Compute distances between all pairs of elements in `itr` according to the
`StringDistance` `dist` and write the result in `r`.

For QGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `itr`. Set `preprocess` to 
false if no preprocessing should be used, regardless of length.
"""
Distances.pairwise!

function Distances.pairwise!(R::AbstractMatrix{<:Number}, dist::StringDistance, X, Y; preprocess = nothing)
    _asymmetric_pairwise!(R, dist, X, Y; preprocess = preprocess)
end

function Distances.pairwise!(R::AbstractMatrix{<:Number}, dist::StringDistance, X; preprocess = nothing)
    (dist isa SemiMetric) ?
        _symmetric_pairwise!(R, dist, X; preprocess = preprocess) :
        _asymmetric_pairwise!(R, dist, X, X; preprocess = preprocess)
end

function _preprocess(X, dist::QGramDistance, preprocess)
    if (preprocess === true) || (isnothing(preprocess) && length(X) >= 5)
        return map(x -> QGramSortedVector(x, dist.q), X)
    else
        return X
    end
end
_preprocess(X, dist::StringDistance, preprocess) = X


function _symmetric_pairwise!(R, dist::StringDistance, X; preprocess = nothing)
    objs = _preprocess(X, dist, preprocess)
    for i in 1:length(objs)
        R[i, i] = 0
        Threads.@threads for j in (i+1):length(objs)
            R[i, j] = R[j, i] = evaluate(dist, objs[i], objs[j])
        end
    end
    return R
end

function _asymmetric_pairwise!(R, dist::StringDistance, X, Y; preprocess = nothing)
    objsX = _preprocess(X, dist, preprocess)
    objsY = _preprocess(Y, dist, preprocess)
    for i in 1:length(objsX)
        Threads.@threads for j in 1:length(objsY)
            R[i, j] = evaluate(dist, objsX[i], objsY[j])
        end
    end
    return R
end