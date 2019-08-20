
using StringDistances, Random
Random.seed!(2)
x = map(Random.randstring, rand(5:25,500_000))
y = map(Random.randstring, rand(5:25,500_000))

function f(t, x, y; min_score = nothing)
    [compare(x[i], y[i], t; min_score = min_score) for i in 1:length(x)]
end

@time f(Hamming(), x, y)
#0.05s
@time f(Jaro(), x, y)
#0.3s
@time f(Levenshtein(), x, y)
# 0.35s. A bit faster than StringDist
@time f(Levenshtein(), x, y, min_score = 0.8)
# 0.11
@time f(DamerauLevenshtein(), x, y)
# 0.45s.  Much faster than StringDist
@time f(DamerauLevenshtein(), x, y, min_score = 0.8)
# 0.08

@time find_best(x[1], y, Levenshtein())
# 0.41
@time find_best(x[1], y, DamerauLevenshtein())
# 0.41

@time find_all(x[1], y, Levenshtein())
# 0.14
@time find_all(x[1], y, DamerauLevenshtein())
# 0.08


# 1.6s slower compared to StringDist






# check
function h(t, x, y; min_score = 1/3)
	out = fill(false, length(x))
	for i in eachindex(x)
		if compare(x[i], y[i], t) <  min_score
			out[i] = compare(x[i], y[i], t ; min_score = min_score) ≈ 0.0
			else
			out[i] = compare(x[i], y[i], t ; min_score = min_score) ≈ compare(x[i], y[i], t)
		end
	end
	all(out)
end
h(Levenshtein(), x, y)
h(DamerauLevenshtein(), x, y)



#= Rcode
library(stringdist)
x <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
system.time(stringdist(x,y,method='lv', nthread = 1))
system.time(stringdist(x,y,method='dl', nthread = 1))

#  0.472
system.time(stringdist(x,y,method='jaccard', nthread = 1))
# 0.739
system.time(stringdist(x,y,method='cosine', nthread = 1))
system.time(stringdist(x,y,method='qgram', nthread = 1))

=#

