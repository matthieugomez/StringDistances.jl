const StringDistance = Union{Hamming, Jaro, JaroWinkler,Levenshtein, OptimalStringAlignement, DamerauLevenshtein, RatcliffObershelp, AbstractQGramDistance, Partial, TokenSort, TokenSet, TokenMax, Normalized}

"""
    compare(s1, s2, dist)

return a similarity score between 0 and 1 for the strings `s1` and 
`s2` based on the distance `dist`.

### Examples
```julia-repl
julia> compare("martha", "marhta", Levenshtein())
0.6666666666666667
```
"""
function compare(s1, s2, dist::StringDistance; min_score = 0.0)
    1 - normalize(dist, max_dist = 1 - min_score)(s1, s2)
end 

"""
    findnearest(s, itr, dist::StringDistance) -> (x, index)

`findnearest` returns the value and index of the element of `itr` that has the 
lowest distance with `s` according to the distance `dist`. 

It is particularly optimized for [`Levenshtein`](@ref) and [`DamerauLevenshtein`](@ref) distances 
(as well as their modifications via [`Partial`](@ref), [`TokenSort`](@ref), [`TokenSet`](@ref), or [`TokenMax`](@ref)).

### Examples
```julia-repl
julia> using StringDistances
julia> s = "Newark"
julia> iter = ["New York", "Princeton", "San Francisco"]
julia> findnearest(s, iter, Levenshtein())
("NewYork", 1)
julia> findnearest(s, iter, Levenshtein(); min_score = 0.9)
(nothing, nothing)
```
"""
function findnearest(s, itr, dist::StringDistance; min_score = 0.0)
    min_score_atomic = Threads.Atomic{Float64}(min_score)
    scores = [0.0 for _ in 1:Threads.nthreads()]
    is = [0 for _ in 1:Threads.nthreads()]
    s = _helper(dist, s)
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, _helper(dist, itr[i]), dist; min_score = min_score_atomic[])
        score_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= score_old
            scores[Threads.threadid()] = score
            is[Threads.threadid()] = i
        end
    end
    imax = is[argmax(scores)]
    imax == 0 ? (nothing, nothing) : (itr[imax], imax)
end
_helper(dist::AbstractQGramDistance, ::Missing) = missing
_helper(dist::AbstractQGramDistance, s) = QGramSortedVector(s, dist.q)
_helper(dist::StringDistance, s) = s

function Base.findmax(s, itr, dist::StringDistance; min_score = 0.0)
    @warn "findmax(s, itr, dist; min_score) is deprecated. Use findnearest(s, itr, dist; min_score)"
    findnearest(s, itr, dist; min_score = min_score)
end

"""
    findall(s, itr , dist::StringDistance; min_score = 0.8)
    
`findall` returns the vector of indices for elements of `itr` that have a 
similarity score higher or equal than `min_score` according to the distance `dist`.
If there are no such elements, return an empty array. 

It is particularly optimized for [`Levenshtein`](@ref) and [`DamerauLevenshtein`](@ref) distances 
(as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`).

### Examples
```julia-repl
julia> using StringDistances
julia> s = "Newark"
julia> iter = ["Newwark", "Princeton", "San Francisco"]
julia> findall(s, iter, Levenshtein())
1-element Array{Int64,1}:
 1
julia> findall(s, iter, Levenshtein(); min_score = 0.9)
0-element Array{Int64,1}
```
"""
function Base.findall(s, itr, dist::StringDistance; min_score = 0.8)
    out = [Int[] for _ in 1:Threads.nthreads()]
    s = _helper(dist, s)
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, _helper(dist, itr[i]), dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end


"""
    pairwise(dist::StringDistance, xs::AbstractVector, ys::AbstractVector = xs; preprocess = nothing)

Compute distances between all pairs of elements in `xs`  and `ys` according to the
`StringDistance` `dist`. Returns a matrix R such that `R[i, j]` corrresponds to the distance between `xs[i]` and `ys[j]`.

For AbstractQGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `xs`. Set `preprocess` to 
false if no preprocessing should be used, regardless of length.

Both symmetric and asymmetric versions are available.

### Examples
```julia-repl
julia> using StringDistances
julia> iter = ["New York", "Princeton"]
julia> pairwise(Levenshtein(), iter)
2×2 Array{Float64,2}:
 0.0  9.0
 9.0  0.0
julia> iter2 = ["San Francisco"]
julia> pairwise(Levenshtein(), iter, iter2)
2×1 Array{Float64,2}:
 12.0
 10.0
```
"""
function pairwise(dist::StringDistance, xs::AbstractVector, ys::AbstractVector = xs; preprocess = nothing)
    T = result_type(dist, eltype(xs), eltype(ys))
    if Missing <: Union{eltype(xs), eltype(ys)}
        T = Union{T, Missing}
    end
    R = Matrix{T}(undef, length(xs), length(ys))
    pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

"""
    pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector, ys::AbstractVector = xs; preprocess = nothing)

Compute distances between all pairs of elements in `xs` and `ys` according to the
`StringDistance` `dist` and write the result in `R`. `R[i, j]` corresponds to the distance between `xs[i]` and `ys[j]`.

For AbstractQGramDistances preprocessing will be used either if `preprocess` is set 
to true or if there are more than 5 elements in `xs`. Set `preprocess` to 
false if no preprocessing should be used, regardless of length.
"""
function pairwise!(R::AbstractMatrix, dist::StringDistance, xs::AbstractVector, ys::AbstractVector = xs; preprocess = nothing)
    length(xs) == size(R, 1) || throw(DimensionMismatch("inconsistent length"))
    length(ys) == size(R, 2) || throw(DimensionMismatch("inconsistent length"))
    ((xs === ys) & (dist isa SemiMetric)) ?
        _symmetric_pairwise!(R, dist, xs; preprocess = preprocess) :
        _asymmetric_pairwise!(R, dist, xs, ys; preprocess = preprocess)
end

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
    objsys = xs === ys ? objsxs : _preprocess(ys, dist, preprocess)
    for i in 1:length(objsxs)
        Threads.@threads for j in 1:length(objsys)
            R[i, j] = evaluate(dist, objsxs[i], objsys[j])
        end
    end
    return R
end

function _preprocess(xs, dist::StringDistance, preprocess)
    if preprocess === nothing
        preprocess = length(xs) >= 5
    end
    if (dist isa AbstractQGramDistance) && preprocess
        return fetch.(map(x -> (Threads.@spawn x === missing ? x : QGramSortedVector(x, dist.q)), xs))
    else
        return xs
    end
end
