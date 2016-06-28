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
- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html)

## Syntax
The function `compare` returns  *a similarity score* between two strings, based on their distance. The similarity score is always between 0 and 1. A value of 0 being completely different and a value of 1 being completely similar.


```julia
using StringDistances
compare(Hamming(), "martha", "marhta")
#> 0.6666666666666667
compare(QGram(2), "martha", "marhta")
#> 0.4
```

To return the *litteral distance* between two strings, use `evaluate`
## Modifiers

The package includes distance modifiers:

- [Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) boosts the similary score of strings with common prefixes

	```julia
	compare(Jaro(), "martha", "marhta")
	#> 0.9444444444444445
	compare(Winkler(Jaro()), "martha", "marhta")
	#> 0.9611111111111111
	```
	The Winkler adjustment was originally defined for the Jaro similarity score but this package defines it for any string distance.

	```julia
	compare(QGram(2), "william", "williams")
	#> 0.9230769230769231
	compare(Winkler(QGram(2)), "william", "williams")
	#> 0.9538461538461539
	```

- The Python library [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) defines a few modifiers for the `RatcliffObershelp` similarity score. This package replicates them and extends them to any string distance:

	- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) returns the maximal similarity score between the shorter string and substrings of the longer string.

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

	- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers by comparing the intersection of two strings with each string.

		```julia
		compare(Jaro(),"mariners vs angels", "los angeles angels at seattle mariners")
		#> 0.559904
		compare(TokenSet(Jaro()),"mariners vs angels", "los angeles angels at seattle mariners")
		#> 0.944444
		```


	- [TokenMax](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) combines scores using the base distance, the `Partial`, `TokenSort` and `TokenSet` modifiers, with penalty terms depending on string lengths.

		```julia
		compare(TokenMax(RatcliffObershelp()),"mariners vs angels", "los angeles angels at seattle mariners")
		#> 0.855
		```

## Unicode
To iterate on graphemes rather than characters, use `graphemeiterator`:

```julia
evaluate(Hamming(), "b\u0300", "a")
#> 2
evaluate(Hamming(), graphemeiterator("b\u0300"), graphemeiterator("a")) 
#> 1
```
## References
- [The stringdist Package for Approximate String Matching](https://journal.r-project.org/archive/2014-1/loo.pdf) Mark P.J. van der Loo
- [fuzzywuzzy blog post](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)


