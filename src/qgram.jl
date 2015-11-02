##############################################################################
##
## Gram Iterator iterates through q-grams of a string
##
## TODO: use Trie? SearchTree?
##############################################################################

type QGramIterator{S <: AbstractString, T <: Integer}
    s::S
    q::T
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
Base.length(qgram::QGramIterator) = max(length(qgram.s) - qgram.q + 1, 0)

function Base.collect(qiter::QGramIterator)
	x = Array(eltype(qiter), length(qiter))
	i = 0
	@inbounds for q in qiter
		i += 1
		x[i] = q
	end
	return x
end

##############################################################################
##
## Define some operations on sorted vector that represent qgrams
##
##############################################################################

function _norm2(v::AbstractVector)
	out = 0
	len = length(v)
	istart = 1
	while istart <= len
		x = v[istart]
		iend = searchsortedlast(v, x, istart, len, Base.Forward)
	    out += (iend - istart + 1)^2
	    istart = iend + 1
	end
	return sqrt(out)
end

function _ndistinct(v::AbstractVector)
	out = 0
	len = length(v)
	istart = 1
	while istart <= len
		x = v[istart]
		iend = searchsortedlast(v, x, istart, len, Base.Forward)
		out += 1
		istart = iend + 1
	end
	return out
end
##############################################################################
##
## q-gram 
## Define v(s) a vector on the space of q-uple which contains number of times it appears in s
## For instance v("leila")["il"] =1 
## q-gram is ∑ |v(s1, p) - v(s2, p)|
##
##############################################################################

type QGram{T <: Integer} <: SemiMetric
    q::T
end
QGram() = QGram(2)

function evaluate(dist::QGram, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer)
    len2 == 0 && return 0

    q1 = QGramIterator(s1, dist.q)
    sort1 = sort!(collect(q1))
    lenq1 = length(sort1)

    q2 = QGramIterator(s2, dist.q)
    sort2 = sort!(collect(q2))
	lenq2 = length(sort2)

	numerator = 0
	i1start = 1
	i2start = 1
	while i1start <= lenq1
		ch1 = sort1[i1start]
		i1end = searchsortedlast(sort1, ch1, i1start, lenq1, Base.Forward)
		i2range = searchsorted(sort2, ch1, i2start, lenq2, Base.Forward)
		numerator += first(i2range) - i2start
		numerator += abs((i1end - i1start + 1) - length(i2range))
		i1start = i1end + 1
		i2start = last(i2range) + 1
	end
	numerator += lenq2 - i2start + 1
    return numerator
end

qgram(s1::AbstractString, s2::AbstractString; q::Integer = 2) = evaluate(QGram(q), s1::AbstractString, s2::AbstractString)

##############################################################################
##
## cosine 
##
## 1 - v(s1, p).v(s2, p)  / ||v(s1, p)|| * ||v(s2, p)||
##############################################################################

type Cosine{T <: Integer} <: SemiMetric
    q::T
end
Cosine() = Cosine(2)


function evaluate(dist::Cosine, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer) 
    len2 == 0 && return 0.0

    q1 = QGramIterator(s1, dist.q)
    sort1 = sort!(collect(q1))
    lenq1 = length(sort1)

    q2 = QGramIterator(s2, dist.q)
	sort2 = sort!(collect(q2))
	lenq2 = length(sort2)

    numerator = 0
    norm1 = 0
    i1start = 1
    i2start = 1
    while i1start <= lenq1
    	ch1 = sort1[i1start]
    	i1end = searchsortedlast(sort1, ch1, i1start, lenq1, Base.Forward)
    	i2range = searchsorted(sort2, ch1, i2start, lenq2, Base.Forward)
        numerator += (i1end - i1start + 1) * length(i2range)
        norm1 += (i1end - i1start + 1)^2
        i1start = i1end + 1
        i2start = last(i2range) + 1
    end

    denominator = sqrt(norm1) * _norm2(sort2)
    return denominator != 0 ? 1.0 - numerator / denominator : s1 == s2 ? 0.0 : 1.0
end

cosine(s1::AbstractString, s2::AbstractString; q::Integer = 2) = evaluate(Cosine(q), s1::AbstractString, s2::AbstractString)

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
type Jaccard{T <: Integer} <: SemiMetric
    q::T
end
Jaccard() = Jaccard(2)

function evaluate(dist::Jaccard, s1::AbstractString, s2::AbstractString, len1::Integer, len2::Integer) 
    len2 == 0 && return 0.0

    q1 = QGramIterator(s1, dist.q)
    sort1 = sort!(collect(q1))
    lenq1 = length(q1)

    q2 = QGramIterator(s2, dist.q)
    sort2 = sort!(collect(q2))
    lenq2 = length(q2)

    numerator = 0
    i1start = 1
    i2start = 1
    norm1 = 0
    while i1start <= lenq1
    	ch1 = sort1[i1start]
    	i1end = searchsortedlast(sort1, ch1, i1start, lenq1, Base.Forward)
    	i2range = searchsorted(sort2, ch1, i2start, lenq2, Base.Forward)
       	numerator += length(i2range) > 0
       	norm1 += 1
       	i1start = i1end + 1
       	i2start = last(i2range) + 1
    end

    norm2 = _ndistinct(sort2)
    denominator = norm1 + norm2 - numerator
    return denominator != 0 ?  1.0 - numerator / denominator : s1 == s2 ? 0.0 : 1.0
end

jaccard(s1::AbstractString, s2::AbstractString; q::Integer = 2) = evaluate(Jaccard(q), s1::AbstractString, s2::AbstractString)