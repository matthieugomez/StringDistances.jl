
using StringDistances, Random
Random.seed!(2)
x = map(Random.randstring, rand(5:25,500_000))
y = map(Random.randstring, rand(5:25,500_000))

function f(t, x, y; max_dist = Inf)
    [evaluate(t, x[i], y[i]; max_dist = max_dist) for i in 1:length(x)]
end

@time f(Hamming(), x, y)
@time f(Jaro(), x, y)
@time f(Levenshtein(), x, y)
# 0.3s. A bit faster than StringDist
@time f(Levenshtein(), x, y, max_dist = 10)
@time f(DamerauLevenshtein(), x, y)
@time f(DamerauLevenshtein(), x, y, max_dist = 10)
# 0.39s.  Much faster than StringDist
@time f(RatcliffObershelp(), x, y)

function g(t, x, y)
    [evaluate(t, x[i], y[i]) for i in 1:length(x)]
end
@time g(Jaccard(2), x, y)
# 1.6s slower compared to StringDist







function h(t, x, y; max_dist = Inf)
    all(evaluate(t, x[i], y[i]; max_dist = max_dist) == min(max_dist, evaluate(t, x[i], y[i])) for i in eachindex(x))
end
h(Jaro(), x, y)
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

