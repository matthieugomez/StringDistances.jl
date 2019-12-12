



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
