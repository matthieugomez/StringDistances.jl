
using DataStructures, StringDistances

x = map(randstring, rand(5:25,100_000))
y = map(randstring, rand(5:25,100_000))
function f(out, t, x, y)
    d = Array(out, length(x))
    @inbounds for i in 1:length(x)
        d[i] = evaluate(t, x[i], y[i])
    end
end

# similar
@time f(Int, Levenshtein(), x, y)
@time f(Float64, Jaro(), x, y)

# 2x slower compared to StringDist
@time f(Int, QGram(2), x, y)
@time f(Float64, Cosine(2), x, y)
@time f(Float64, Jaccard(2), x, y)

#
@time f(Float64, RatcliffObershelp(), x, y)






#= Rcode
library(stringdist)
x <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
system.time(stringdist(x,y,method='lv', nthread = 1))
system.time(stringdist(x,y,method='jaccard', nthread = 1))
system.time(stringdist(x,y,method='cosine', nthread = 1))
system.time(stringdist(x,y,method='qgram', nthread = 1))

=#

