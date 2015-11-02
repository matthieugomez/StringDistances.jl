
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

# all 5-10x slower compared to StringDist
@time f(Int, QGram(2), x, y)
@time f(Float64, Cosine(2), x, y)
@time f(Float64, Jaccard(2), x, y)






#= Rcode
library(stringdist)
x <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse="")) 
y <- sapply(sample(5:25,1e5,replace=TRUE), function(n) paste(sample(letters,n,replace=TRUE),collapse=""))
system.time(stringdist(x,y,method='lv'))
system.time(stringdist(x,y,method='jaccard'))
system.time(stringdist(x,y,method='cosine'))
system.time(stringdist(x,y,method='qgram'))

=#



function f(x, y)
    d = Array(Float64, length(x))
    sort1 = Array(SubString{ASCIIString}, 25)
    sort2 = Array(SubString{ASCIIString}, 25)
    @inbounds for i in 1:length(x)
        d[i] = evaluate(Jaccard(2), x[i], y[i], sort1 sort2)
    end
end
@time f(x, y)




function g(x, y)
    d = Array(Float64, length(x))
    set1 = Set{SubString{ASCIIString}}()
    set2 = Set{SubString{ASCIIString}}()
    @inbounds for i in 1:length(x)
        d[i] = evaluate(Jaccard(2), x[i], y[i], set1, set2)
    end
end
@time g(x, y)

