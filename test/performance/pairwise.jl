using StringDistances, Random
using BenchmarkTools

N = if length(ARGS) > 0
    try
        parse(Int, ARGS[1])
    catch _
        100
    end
else
    100 # default value
end

Maxlength = if length(ARGS) > 1
    try
        parse(Int, ARGS[2])
    catch _
        100
    end
else
    100 # default value
end

S = String[randstring(rand(3:Maxlength)) for _ in 1:N]

println("For ", Threads.nthreads(), " threads and ", N, " strings of max length ", Maxlength, ":")

dist = Cosine(2)
t1 = @belapsed dm1 = pairwise(dist, S; preprocess = false)
t2 = @belapsed dm2 = pairwise(dist, S; preprocess = true)

println("  - time WITHOUT pre-calculation: ", round(t1, digits = 3))
println("  - time WITH    pre-calculation: ", round(t2, digits = 3))
println("  - speedup with pre-calculation: ", round(t1/t2, digits = 1))
