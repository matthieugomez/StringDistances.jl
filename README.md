[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)
[![StringDistances](http://pkg.julialang.org/badges/StringDistances_0.4.svg)](http://pkg.julialang.org/?pkg=StringDistances)

This Julia package computes various distances between strings.



## Distances

#### Edit Distances
- Hamming Distance
- Jaro Distance
- Levenshtein Distance
- Damerau-Levenshtein Distance
- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html) (similar to the Python library [difflib](https://docs.python.org/2/library/difflib.html))

#### Q-Grams Distances
- QGram Distance
- Cosine Distance
- Jaccard Distance

A good reference for q-gram distances is the article written for the R package `stringdist`:
*The stringdist Package for Approximate String Matching* Mark P.J. van der Loo


## Syntax



#### evaluate
The function `evaluate` returns the litteral distance between two strings (a value of 0 being identical). While some distances are bounded by 1, other distances like `Hamming`, `Levenshtein`, `Damerau-Levenshtein`,  `Jaccard` can be higher than 1.

```julia
using StringDistances
evaluate(Hamming(), "martha", "marhta")
#> 2
evaluate(QGram(2), "martha", "marhta")
#> 6
```

#### compare
The higher level function `compare` directly computes for any distance a similarity score between 0 and 1. A value of 0 being completely different and a value of 1 being completely similar.
```julia
using StringDistances
compare(Hamming(), "martha", "marhta")
#> 0.6666666666666667
compare(QGram(2), "martha", "marhta")
#> 0.4
```


## Modifiers

The package defines a number of types to modify string metrics:

- [Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) boosts the similary score of strings with common prefixes

	```julia
	compare(Jaro(), "martha", "marhta")
	#> 0.9444444444444445
	compare(Winkler(Jaro()), "martha", "marhta")
	#> 0.9611111111111111
	```
	The Winkler adjustment was originally defined for the Jaro distance but this package defines it for any string distance.

	```julia
	compare(QGram(2), "william", "williams")
	#> 0.9230769230769231
	compare(Winkler(QGram(2)), "william", "williams")
	#> 0.9538461538461539
	```

- For strings composed of several words, the Python library [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) defines a few modifiers for the `RatcliffObershelp` distance. This package defines them for any string distance:

	- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in string lengths. The function returns the maximal similarity score between the shorter string and all substrings of the longer string. 	

		```julia
		compare(Partial(Hamming()), "New York Yankees", "Yankees")
		#> 1.0
		```

	- [TokenSort](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders by reording words alphabetically.

		```julia
		compare(TokenSort(RatcliffObershelp()),"mariners vs angels", "angels vs mariners")
		#> 1.0
		```

	- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers.

		```julia
		compare(TokenSet(RatcliffObershelp()),"mariners vs angels", "los angeles angels of anaheim at seattle mariners")
		```


