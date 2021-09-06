
using StringDistances, Random
Random.seed!(2)
x = map(Random.randstring, rand(5:25,500_000))
y = map(Random.randstring, rand(5:25,500_000))

function f(t, x, y; min_score = 0.0)
    [compare(x[i], y[i], t; min_score = min_score) for i in 1:length(x)]
end

function g(dist, x, y)
    [dist(x[i], y[i]) for i in 1:length(x)]
end



@time f(Jaro(), x, y);
#0.3s 
@time f(Levenshtein(), x, y);
# 0.4s
@time f(Levenshtein(), x, y, min_score = 0.8);
# 0.11 
@time f(DamerauLevenshtein(), x, y);
# 0.58s.
@time f(DamerauLevenshtein(), x, y, min_score = 0.8);
# 0.08 (now 0.09)
@time f(RatcliffObershelp(), x, y);
# 1.35s




@time findnearest(x[1], y, Levenshtein());
# 0.1
@time findnearest(x[1], y, DamerauLevenshtein());
# 0.1
@time findnearest(x[1], y, QGram(2));
# 0.75



@time findall(x[1], y, Levenshtein());
# 0.05
@time findall(x[1], y, DamerauLevenshtein());
# 0.05
@time findall(x[1], y, Partial(DamerauLevenshtein()));
# 0.96
@time findall(x[1], y, QGram(2));
# 0.81
@time findall(x[1], y, TokenSort(DamerauLevenshtein()));
# 0.27 (now 0.32)
@time findall(x[1], y, TokenSet(DamerauLevenshtein()));
# 0.55
@time findall(x[1], y, TokenMax(DamerauLevenshtein()));
# 2.25 (now 3.6)


x = map(Random.randstring, rand(5:25,1000))
y = map(Random.randstring, rand(5:25,1000))
@time pairwise(Levenshtein(), x, y)
# 0.25 seconds
@time pairwise(QGram(2), x, y, preprocess = false)
# 2.126829 
@time pairwise(QGram(2), x, y, preprocess = true)
# 0.12




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

