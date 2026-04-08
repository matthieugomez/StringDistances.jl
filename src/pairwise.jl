"""
    pairwise(dist::StringDistance, xs::AbstractVector, ys::AbstractVector = xs; preprocess = true)

Compute distances between all pairs of elements in `xs`  and `ys` according to the
`StringDistance` `dist`. Returns a matrix R such that `R[i, j]` corrresponds to the distance between `xs[i]` and `ys[j]`.

Set `preprocess` to false if no preprocessing should be used.

Both symmetric and asymmetric versions are available.

### Examples
```julia-repl
julia> using StringDistances
julia> iter = ["New York", "Princeton"]
julia> pairwise(Levenshtein(), iter)
2×2 Matrix{Int64}:
 0  9
 9  0
julia> iter2 = ["San Francisco"]
julia> pairwise(Levenshtein(), iter, iter2)
2×1 Matrix{Int64}:
 12
 10
```
"""
function StatsAPI.pairwise(dist::Union{StringSemiMetric, StringMetric}, xs::AbstractVector, ys::AbstractVector = xs; preprocess = true)
    T = result_type(dist, eltype(xs), eltype(ys))
    R = Matrix{T}(undef, length(xs), length(ys))
    pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

"""
    pairwise!(R::AbstractMatrix, dist::Union{StringSemiMetric, StringMetric}, xs::AbstractVector, ys::AbstractVector = xs; preprocess = true)

Compute distances between all pairs of elements in `xs` and `ys` according to the
`Union{StringSemiMetric, StringMetric}` `dist` and write the result in `R`. `R[i, j]` corresponds to the distance between `xs[i]` and `ys[j]`.

Set `preprocess` to false if no preprocessing should be used.
"""
function StatsAPI.pairwise!(R::AbstractMatrix, dist::Union{StringSemiMetric, StringMetric}, xs::AbstractVector, ys::AbstractVector = xs; preprocess = true)
    length(xs) == size(R, 1) || throw(DimensionMismatch("inconsistent length"))
    length(ys) == size(R, 2) || throw(DimensionMismatch("inconsistent length"))
    (xs === ys) ?
        _symmetric_pairwise!(R, dist, xs; preprocess = preprocess) :
        _asymmetric_pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

function _symmetric_pairwise!(R::AbstractMatrix, dist::Union{StringSemiMetric, StringMetric}, xs::AbstractVector; preprocess = true)
    if preprocess
        xs = _preprocess_list(dist, xs)
    end
    Threads.@threads :dynamic for i in 1:length(xs)
        # handle missing
        R[i, i] = xs[i] != xs[i]
        for j in (i+1):length(xs)
            R[i, j] = R[j, i] = evaluate(dist, xs[i], xs[j])
        end
    end
    return R
end

function _asymmetric_pairwise!(R::AbstractMatrix, dist::Union{StringSemiMetric, StringMetric}, xs::AbstractVector, ys::AbstractVector; preprocess = true)
    if preprocess
        objxs = _preprocess_list(dist, xs)
        objys = xs === ys ? objxs : _preprocess_list(dist, ys)
    else
        objxs = xs
        objys = ys
    end
    if length(objxs) >= length(objys)
        Threads.@threads :dynamic for i in 1:length(objxs)
            for j in 1:length(objys)
                R[i, j] = evaluate(dist, objxs[i], objys[j])
            end
        end
    else
        Threads.@threads :dynamic for j in 1:length(objys)
            for i in 1:length(objxs)
                R[i, j] = evaluate(dist, objxs[i], objys[j])
            end
        end
    end
    return R
end

_preprocess_list(dist::Union{StringSemiMetric, StringMetric}, xs)  = xs

function _preprocess_list(dist::AbstractQGramDistance, xs::AbstractVector)
    idx = findfirst(x -> x !== missing, xs)
    idx === nothing && return copy(xs)

    S = Base.nonmissingtype(eltype(xs))
    if !isconcretetype(S)
        return map(x -> x === missing ? x : QGramSortedVector(x, dist.q), xs)
    end

    sample = QGramSortedVector(xs[idx], dist.q)
    T = Missing <: eltype(xs) ? Union{Missing, typeof(sample)} : typeof(sample)
    out = Vector{T}(undef, length(xs))
    ranges = _chunk_ranges(length(xs))

    Threads.@threads :dynamic for ir in eachindex(ranges)
        for i in ranges[ir]
            x = xs[i]
            out[i] = x === missing ? missing : QGramSortedVector(x, dist.q)
        end
    end
    return out
end
