@doc """
    pairwise(dist::StringDistance, xs::AbstractVector; preprocess = nothing)
    pairwise(dist::StringDistance, xs::AbstractVector, ys::AbstractVector; preprocess = nothing)

Compute distances between all pairs of elements in `xs`  and `ys` according to the
`StringDistance` `dist`.

For QGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `xs`. Set `preprocess` to 
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

function Distances.pairwise(dist::StringDistance, xs::AbstractVector, ys::AbstractVector; preprocess = nothing)
    T = result_type(dist, eltype(xs), eltype(ys))
    if Missing <: Union{eltype(xs), eltype(ys)}
        T = Union{T, Missing}
    end
    R = Matrix{T}(undef, length(xs), length(ys))
    pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

function Distances.pairwise(dist::StringDistance, xs::AbstractVector; preprocess = nothing)
    T = result_type(dist, eltype(xs), eltype(xs))
    if Missing <: eltype(xs)
        T = Union{T, Missing}
    end
    R = Matrix{T}(undef, length(xs), length(xs))
    pairwise!(R, dist, xs; preprocess = preprocess)
end

@doc """
    pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector; preprocess = nothing)
    pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector, ys::AbstractVector; preprocess = nothing)

Compute distances between all pairs of elements in `xs` and `ys` according to the
`StringDistance` `dist` and write the result in `R`.

For QGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `xs`. Set `preprocess` to 
false if no preprocessing should be used, regardless of length.
"""
Distances.pairwise!

function Distances.pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector, ys::AbstractVector; preprocess = nothing)
    length(xs) == size(R, 1) || throw(DimensionMismatch("inconsistent length"))
    length(ys) == size(R, 2) || throw(DimensionMismatch("inconsistent length"))
    _asymmetric_pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

function Distances.pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector; preprocess = nothing)
    length(xs) == size(R, 1) || throw(DimensionMismatch("inconsistent length"))
    length(xs) == size(R, 2) || throw(DimensionMismatch("inconsistent length"))
    (dist isa SemiMetric) ?
        _symmetric_pairwise!(R, dist, xs; preprocess = preprocess) :
        _asymmetric_pairwise!(R, dist, xs, xs; preprocess = preprocess)
end

function _preprocess(xs, dist::QGramDistance, preprocess)
    if preprocess === nothing ? length(xs) >= 5 : preprocess 
        return map(x -> x === missing ? x : QGramSortedVector(x, dist.q), xs)
    else
        return xs
    end
end
_preprocess(xs, dist::StringDistance, preprocess) = xs


function _symmetric_pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector; preprocess = nothing)
    objs = _preprocess(xs, dist, preprocess)
    for i in 1:length(objs)
        # handle missing
        R[i, i] = objs[i] != objs[i]
        Threads.@threads for j in (i+1):length(objs)
            R[i, j] = R[j, i] = evaluate(dist, objs[i], objs[j])
        end
    end
    return R
end

function _asymmetric_pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector, ys::AbstractVector; preprocess = nothing)
    objsxs = _preprocess(xs, dist, preprocess)
    objsys = _preprocess(ys, dist, preprocess)
    for i in 1:length(objsxs)
        Threads.@threads for j in 1:length(objsys)
            R[i, j] = evaluate(dist, objsxs[i], objsys[j])
        end
    end
    return R
end
