
using StringDistances
srand(2)
x = map(Base.randstring, rand(5:25,500_000))
y = map(Base.randstring, rand(5:25,500_000))
function f(t, x, y)
    [evaluate(t, x[i], y[i]) for i in 1:length(x)]
end

# same speed as StringDist
@time f(Levenshtein(), x, y)
@time f(Jaro(), x, y)
@time f(RatcliffObershelp(), x, y)

# 4x slower compared to StringDist
@time f(Jaccard(2), x, y)
@time f(Cosine(2), x, y)
@time f(QGram(2), x, y)

#






#= Rcode
library(stringdist)
x <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,5 * 1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
system.time(stringdist(x,y,method='lv', nthread = 1))
system.time(stringdist(x,y,method='jaccard', nthread = 1))
system.time(stringdist(x,y,method='cosine', nthread = 1))
system.time(stringdist(x,y,method='qgram', nthread = 1))

=#

