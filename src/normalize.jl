struct Normalized{V <: SemiMetric} <: SemiMetric
    dist::V
    max_dist::Float64
end

function (dist::Normalized{<:Hamming})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    out = dist.dist(s1, s2) / len2
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:Union{Levenshtein{Nothing}, DamerauLevenshtein{Nothing}}})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len2 == 0 && return 1.0
    if dist.dist isa Levenshtein
        d = Levenshtein(ceil(Int, len2 * dist.max_dist))(s1, s2)
    else
        d = DamerauLevenshtein(ceil(Int, len2 * dist.max_dist))(s1, s2)
    end
    out = d / len2
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized{<:AbstractQGramDistance})(s1, s2)
    ((s1 === missing) | (s2 === missing)) && return missing
    # When string length < q for qgram distance, returns s1 == s2
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    len1 <= dist.dist.q - 1 && return convert(Float64, s1 != s2)
    if dist.dist isa QGram
        out = dist.dist(s1, s2) / (len1 + len2 - 2 * dist.dist.q + 2)
    else
        out = dist.dist(s1, s2)
    end
    out > dist.max_dist ? 1.0 : out
end

function (dist::Normalized)(s1, s2)
    out = dist.dist(s1, s2)
    out > dist.max_dist ? 1.0 : out
end

"""
   normalize(dist)

Creates a normalized distance. The distance always return a Float64 between 0.0 and 1.0 (or a missing if one of the argument is missing)

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> Levenshtein()(s1, s2)
25
julia> StringDistances.normalize(Levenshtein())(s1, s2)
0.8064 
```
"""
normalize(dist::SemiMetric; max_dist = 1.0) = Normalized{typeof(dist)}(dist, max_dist)
normalize(dist::Partial; max_dist = 1.0) = Partial(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSort; max_dist = 1.0) = TokenSort(normalize(dist.dist; max_dist = max_dist))
normalize(dist::TokenSet; max_dist = 1.0) = TokenSet(normalize(dist.dist; max_dist = max_dist))
normalize(dist::Normalized; max_dist = 1.0) = Normalized(dist.dist, max_dist)


"""
   TokenMax(dist)

Creates the `TokenMax{dist}` distance

`TokenMax{dist}` normalizes the distance `dist` and returns the minimum of the distance,
its [`Partial`](@ref) modifier, its [`TokenSort`](@ref) modifier, and its 
[`TokenSet`](@ref) modifier, with penalty terms depending on string length.

It is only defined on AbstractStrings

### Examples
```julia-repl
julia> s1 = "New York Mets vs Atlanta"
julia> s2 = "Atlanta Braves vs New York Mets"
julia> evaluate(TokenMax(RatcliffObershelp()), s1, s2)
0.05
```
"""
struct TokenMax{S <: SemiMetric} <: SemiMetric
    dist::S
    max_dist::Float64
end
TokenMax(dist::SemiMetric; max_dist = 1.0) = TokenMax(dist, max_dist)
normalize(dist::TokenMax; max_dist = 1.0) = TokenMax(dist.dist, max_dist)

function (dist::TokenMax)(s1::Union{AbstractString, Missing}, s2::Union{AbstractString, Missing})
    ((s1 === missing) | (s2 === missing)) && return missing
    s1, s2 = reorder(s1, s2)
    len1, len2 = length(s1), length(s2)
    max_dist = dist.max_dist
    dist0 = normalize(dist.dist; max_dist = max_dist)
    score = dist0(s1, s2)
    min_score = min(max_dist, score)
    unbase_scale = 0.95
    # if one string is much shorter than the other, use partial
    if length(s2) >= 1.5 * length(s1)
        partial_scale = length(s2) > (8 * length(s1)) ? 0.6 : 0.9
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / partial_scale)
        score_partial = 1 - partial_scale * (1 - Partial(dist0)(s1, s2))
        min_score = min(max_dist, score_partial)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_sort = 1 - unbase_scale * partial_scale * (1 - TokenSort(Partial(dist0))(s1, s2))
        max_dist = min(max_dist, score_sort)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / (unbase_scale * partial_scale))
        score_set = 1 - unbase_scale * partial_scale * (1 - TokenSet(Partial(dist0))(s1, s2)) 
        out = min(score, score_partial, score_sort, score_set)
    else
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / unbase_scale)
        score_sort = 1 - unbase_scale * (1 - TokenSort(dist0)(s1, s2))
        max_dist = min(max_dist, score_sort)
        dist0 = normalize(dist0, max_dist = 1 - (1 - max_dist) / unbase_scale)
        score_set = 1 - unbase_scale * (1 - TokenSet(dist0)(s1, s2))
        out = min(score, score_sort, score_set)
    end
    out > max_dist ? 1.0 : out
end




const StringDistance = Union{Hamming, Jaro, JaroWinkler,Levenshtein, DamerauLevenshtein, RatcliffObershelp, AbstractQGramDistance, Partial, TokenSort, TokenSet, TokenMax, Normalized}

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
    _findnearest(s, itr, dist; min_score = min_score)
end
function findnearest(s, itr, dist::AbstractQGramDistance; min_score = 0.0)
    _findnearest(QGramSortedVector(s, dist.q), itr, dist; min_score = min_score)
end

function _findnearest(s, itr, dist::StringDistance; min_score = 0.0)
    min_score_atomic = Threads.Atomic{Float64}(min_score)
    scores = [0.0 for _ in 1:Threads.nthreads()]
    is = [0 for _ in 1:Threads.nthreads()]
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, itr[i], dist; min_score = min_score_atomic[])
        score_old = Threads.atomic_max!(min_score_atomic, score)
        if score >= score_old
            scores[Threads.threadid()] = score
            is[Threads.threadid()] = i
        end
    end
    imax = is[argmax(scores)]
    imax == 0 ? (nothing, nothing) : (itr[imax], imax)
end


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
    _findall(s, itr, dist; min_score = min_score)
end
function Base.findall(s, itr, dist::AbstractQGramDistance; min_score = 0.8)
    _findall(QGramSortedVector(s, dist.q), itr, dist; min_score = min_score)
end
function _findall(s, itr, dist::StringDistance; min_score = 0.8)
    out = [Int[] for _ in 1:Threads.nthreads()]
    # need collect since @threads requires a length method
    Threads.@threads for i in collect(eachindex(itr))
        score = compare(s, itr[i], dist; min_score = min_score)
        if score >= min_score
            push!(out[Threads.threadid()], i)
        end
    end
    vcat(out...)
end
