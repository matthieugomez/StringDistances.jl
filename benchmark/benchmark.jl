
using StringDistances

x = map(randstring, rand(5:25,100_000))
y = map(randstring, rand(5:25,100_000))
function f(out, t, x, y)
    d = Array(out, length(x))
    for i in 1:length(x)
        d[i] = StringDistances.evaluate(t, x[i], y[i])
    end
end

# I get 0.12 vs 0.10 in stringdist
# http://www.markvanderloo.eu/yaRb/2013/09/07/a-bit-of-benchmarking-with-string-distances/

@time f(Int, Levenshtein(), x, y)
@time f(Float64, Jaccard(2), x, y)
@time f(Float64, Cosine(2), x, y)



#= Rcode
library(stringdist)
x <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
stringdist(x,y,method='lv')
stringdist(x,y,method='jaccard')
stringdist(x,y,method='jaccard')
stringdist(x,y,method='cosine')

=#