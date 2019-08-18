
using StringDistances, Random
Random.seed!(2)
x = map(Random.randstring, rand(5:25,500_000))
y = map(Random.randstring, rand(5:25,500_000))
function f(t, x, y)
    [evaluate(t, x[i], y[i]) for i in 1:length(x)]
end
@time f(Hamming(), x, y)
@time f(Jaro(), x, y)
@time f(Levenshtein(), x, y)
# 0.3s. A big faster than StringDist
@time f(DamerauLevenshtein(), x, y)
@time f(RatcliffObershelp(), x, y)
@time f(Jaccard(2), x, y)
# 1.6s 2-3x slower compared to StringDist

# a bist faster than StringDist
@time f(Levenshtein(), x, y)
#  355.984 ms (1500004 allocations: 223.24 MiB)
@time f(RatcliffObershelp(), x, y)

@time f(Jaccard(2), x, y)
# 1.6s






#= Rcode
library(stringdist)
x <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
system.time(stringdist(x,y,method='lv', nthread = 1))
#  0.472
system.time(stringdist(x,y,method='jaccard', nthread = 1))
# 0.739
system.time(stringdist(x,y,method='cosine', nthread = 1))
system.time(stringdist(x,y,method='qgram', nthread = 1))

=#

