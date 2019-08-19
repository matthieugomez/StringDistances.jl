function h(t, x, y; max_dist = Inf)
    all(evaluate(t, x[i], y[i]; max_dist = max_dist) == min(max_dist, evaluate(t, x[i], y[i])) for i in eachindex(x))
end
h(Jaro(), x, y)
h(Levenshtein(), x, y)
h(DamerauLevenshtein(), x, y)