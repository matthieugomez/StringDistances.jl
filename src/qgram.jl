##############################################################################
##
## Define QGram Distance type
##
##############################################################################
abstract AbstractQGram <: SemiMetric

##############################################################################
##
## Define a type that iterates through q-grams of a string
##
##############################################################################

type QGramIterator{S <: AbstractString, T <: Integer}
	s::S # string
	l::Int # length of string
	q::T # length of q-grams
end
function Base.start(qgram::QGramIterator)
	len = length(qgram.s)
	(1, len < qgram.q ? endof(qgram.s) + 1 : chr2ind(qgram.s, qgram.q))
end
function Base.next(qgram::QGramIterator, state)
	istart, iend = state
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	return element, nextstate
end
function Base.done(qgram::QGramIterator, state)
	istart, idend = state
	done(qgram.s, idend)
end
Base.eltype(qgram::QGramIterator) = SubString{typeof(qgram.s)}
Base.length(qgram::QGramIterator) = max(qgram.l - qgram.q + 1, 0)
function Base.collect(qgram::QGramIterator)
	x = Array(eltype(qgram), length(qgram))
	i = 0
	for q in qgram
		i += 1
		@inbounds x[i] = q
	end
	return x
end
Base.sort(qgram::QGramIterator) = sort!(collect(qgram), alg = QuickSort)

##############################################################################
##
## Define a Tree
##
##############################################################################

abstract Tree{K, V}
type EmptyTree{K,V} <: Tree{K, V}
end

type TreeNode{K,V} <: Tree{K, V}
    key::  K
    data::Tuple{V, V}
    left:: Tree{K,V}
    right::Tree{K,V}
end

add1!{K,V}(t::EmptyTree{K,V}, k) = TreeNode{K,V}(k, (one(V), zero(V)), t, t)
function add1!{K, V}(t::TreeNode{K, V}, k)
    if t.key == k
    	a, b = t.data
        t.data = (a + one(V), b)
    elseif k < t.key
        t.left = add1!(t.left, k)
    else
        t.right = add1!(t.right, k)
    end
    return t
end

add2!{K,V}(t::EmptyTree{K,V}, k) = TreeNode{K,V}(k, (zero(V), one(V)), t, t)
function add2!{K, V}(t::TreeNode{K, V}, k)
    if t.key == k
        a, b = t.data
        t.data = (a, b + one(V))
    elseif k < t.key
        t.left = add2!(t.left, k)
    else
        t.right = add2!(t.right, k)
    end
    return t
end

function Tree{S}(dist, s1::S, s2::S, len1::Integer, len2::Integer)
	qgram1 = QGramIterator(s1, len1, dist.q)
	qgram2 = QGramIterator(s2, len2, dist.q)
	t = EmptyTree{SubString{S}, UInt}()
	for x in qgram1
		t = add1!(t, x)
	end
	for x in qgram2
		t = add2!(t, x) 
	end
	return t
end


##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

type QGram{T <: Integer} <: AbstractQGram
	q::T
end
QGram() = QGram(2)

type QGramAccumulator
	n::Int
end

function evaluate(dist::QGram, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
	isempty(s2) && return 0
	t = Tree(dist, s1, s2, len1, len2)
	acc = QGramAccumulator(0)
	evalans(dist, t, acc)
	return acc.n
end
evalans(dist::QGram, t::EmptyTree, acc::QGramAccumulator) = nothing
function evalans{K, V}(dist::QGram, t::TreeNode{K, V}, acc::QGramAccumulator)
	n1, n2 = t.data
	acc.n += n1 > n2 ? n1 - n2 : n2 - n1
	evalans(dist, t.left, acc)
	evalans(dist, t.right, acc)
end

##############################################################################
##
## cosine 
##
## 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
##############################################################################

type Cosine{T <: Integer} <: AbstractQGram
	q::T
end
Cosine() = Cosine(2)

type CosineAccumulator
	norm1::Int
	norm2::Int
	prodnorm::Int
end
function evaluate(dist::Cosine, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
	isempty(s2) && return 0.0
	t = Tree(dist, s1, s2, len1, len2)
	acc = CosineAccumulator(0, 0, 0)
	evalans(dist, t, acc)	
	denominator = sqrt(acc.norm1) * sqrt(acc.norm2)
	return denominator != 0 ? 1.0 - acc.prodnorm / denominator : s1 == s2 ? 0.0 : 1.0
end

evalans(dist::Cosine, t::EmptyTree, acc::CosineAccumulator) = nothing
function evalans{K, V}(dist::Cosine, t::TreeNode{K, V}, acc::CosineAccumulator)
	n1, n2 = t.data
	acc.norm1 += n1^2
	acc.norm2 += n2^2
	acc.prodnorm += n1 * n2
	evalans(dist, t.left, acc)
	evalans(dist, t.right, acc)
end
function cosine(s1::AbstractString, s2::AbstractString; q::Integer = 2)
	evaluate(Cosine(q), s1::AbstractString, s2::AbstractString)
end

##############################################################################
##
## Jaccard
##
## Denote Q(s, q) the set of tuple of length q in s
## 1 - |intersect(Q(s1, q), Q(s2, q))| / |union(Q(s1, q), Q(s2, q))|
##
## return 1.0 if smaller than qgram
##
##############################################################################
type Jaccard{T <: Integer} <: AbstractQGram
	q::T
end
Jaccard() = Jaccard(2)

type JaccardAccumulator
	ndistinct1::Int
	ndistinct2::Int
	nintersect::Int
end

function evaluate(dist::Jaccard, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
	isempty(s2) && return 0.0
	t = Tree(dist, s1, s2, len1, len2)
	acc = JaccardAccumulator(0, 0, 0)
	evalans(dist, t, acc)
	denominator = acc.ndistinct1 + acc.ndistinct2 - acc.nintersect
	return denominator != 0 ? 1.0 - acc.nintersect / denominator : s1 == s2 ? 0.0 : 1.0
end


evalans(dist::Jaccard, t::EmptyTree, acc::JaccardAccumulator) = nothing
function evalans{K, V}(dist::Jaccard, t::TreeNode{K, V}, acc::JaccardAccumulator)
	n1, n2 = t.data
	acc.ndistinct1 += (n1 > zero(V))
	acc.ndistinct2 += (n2 > zero(V))
	acc.nintersect += ((n1 > zero(V)) & (n2 > zero(V)))
	evalans(dist, t.left, acc)
	evalans(dist, t.right, acc)
end

function jaccard(s1::AbstractString, s2::AbstractString; q::Integer = 2)
	evaluate(Jaccard(q), s1::AbstractString, s2::AbstractString)
end


