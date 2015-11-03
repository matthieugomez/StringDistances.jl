[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)
[![StringDistances](http://pkg.julialang.org/badges/StringDistances_0.4.svg)](http://pkg.julialang.org/?pkg=StringDistances)

StringDistances allow to compute various distances between strings. The package should work with any `AbstractString` (in particular ASCII and UTF-8)


## Distances

- Hamming Distance
- Jaro Distance
- Levenshtein Distance
- Damerau-Levenshtein Distance
- QGram Distance
- Cosine Distance
- Jaccard Distance


A good reference about string distances is the article written for the R package `stringdist`:
*The stringdist Package for Approximate String Matching* Mark P.J. van der Loo

## Syntax
- The basic syntax follows the [Distances](https://github.com/JuliaStats/Distances.jl) package:

	```julia
	using StringDistances
	evaluate(Hamming(), "martha", "marhta")
	evaluate(QGram(2), "martha", "marhta")
	```

- Normalize a distance between 0-1 with `Normalized`

	```julia
	evaluate(Normalized(Hamming()), "martha", "marhta")
	evaluate(Normalized(QGram(2)), "martha", "marhta")
	```

- Add a [Winkler adjustment](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) with `Winkler`

	```julia
	evaluate(Winkler(Jaro()), "martha", "marhta")
	evaluate(Winkler(Qgram(2)), "martha", "marhta")
	```
	While the Winkler adjustment was originally defined in the context of the Jaro distance, it can be helpful with other distances too. Note: a distance is automatically normalized between 0 and 1 when used with a Winkler adjustment.
