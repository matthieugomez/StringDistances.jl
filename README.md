[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)
[![StringDistances](http://pkg.julialang.org/badges/StringDistances_0.4.svg)](http://pkg.julialang.org/?pkg=StringDistances)

This Julia package computes various distances between strings.



## Distances

#### Edit Distances
- [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance)
- [Levenshtein Distance](https://en.wikipedia.org/wiki/Levenshtein_distance)
- [Damerau-Levenshtein Distance](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance)

#### Q-Grams Distances
Q-gram distances compare the set of all substrings of length `q` in each string.
- QGram Distance
- [Cosine Distance](https://en.wikipedia.org/wiki/Cosine_similarity)
- [Jaccard Distance](https://en.wikipedia.org/wiki/Jaccard_index)
- [Overlap Distance](https://en.wikipedia.org/wiki/Overlap_coefficient)
- [Sorensen-Dice Distance](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient)

#### Others
- [Jaro Distance](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance)
- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html) is based on the length of matching subsequences. It is used in the Python library [difflib](https://docs.python.org/2/library/difflib.html).

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

The package defines a number of ways to modify string metrics:

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

- The Python library [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) defines a few modifiers for the `RatcliffObershelp` distance. This package defines them for any string distance:

	- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in string lengths. The function returns the maximal similarity score between the shorter string and all substrings of the longer string. 	

		```julia
		compare(Levenshtein(), "New York Yankees", "Yankees")
		#> 0.4375
		compare(Partial(Levenshtein()), "New York Yankees", "Yankees")
		#> 1.0
		```

	- [TokenSort](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders by reording words alphabetically.

		```julia
		compare(RatcliffObershelp(), "mariners vs angels", "angels vs mariners")
		#> 0.44444
		compare(TokenSort(RatcliffObershelp()),"mariners vs angels", "angels vs mariners")
		#> 1.0
		```

	- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers.

		```julia
		compare(Jaro(),"mariners vs angels", "los angeles angels at seattle mariners")
		#> 0.559904
		compare(TokenSet(Jaro()),"mariners vs angels", "los angeles angels at seattle mariners")
		#> 0.944444
		```


- You can compose multiple modifiers:
	```julia
	compare(Winkler(Partial(Jaro())),"mariners vs angels", "los angeles angels at seattle mariners")
	#> 0.7378917378917379
	compare(TokenSet(Partial(RatcliffObershel())),"mariners vs angels", "los angeles angels at seattle mariners")
	#> 1.0
	```


## Tips
In case you're wondering which distance to use:

- Each distance is tailored to a specific problem. Edit distances works well with local spelling errors, the Ratcliff-Obsershelp distance works well with edited texts, the Jaro Winkler distance was invented for short strings such as person names, the QGrams distances works well with strings composed of multiple words with fluctuating orderings.
- When comparing company or individual names, each string is composed of multiple words and their ordering is mostly irrelevant. Edit distances will perform poorly in this situation. Use either a distance robust to word order (like QGram distances), or compose a distance with `TokenSort` or `TokenSet`, which reorder the words alphabetically.

	```julia
	compare(RatcliffObershelp(), "mariners vs angels", "angels vs mariners")
	#> 0.44444
	compare(TokenSort(RatcliffObershelp()),"mariners vs angels", "angels vs mariners")
	#> 1.0
	compare(Cosine(3), "mariners vs angels", "angels vs mariners")
	#> 0.8125
	```

- Standardize strings before comparing them (lowercase, punctuation, whitespaces, accents, abbreviations...)



## References
- [The stringdist Package for Approximate String Matching](https://journal.r-project.org/archive/2014-1/loo.pdf) Mark P.J. van der Loo
- [fuzzywuzzy blog post](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)


