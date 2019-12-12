[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)

This Julia package computes various distances between `AbstractString`s

## Installation
The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and so can be installed at the REPL with `] add StringDistances`.

## Syntax
The function `compare` returns a similarity score between two strings. The function always returns a score between 0 and 1, with a value of 0 being completely different and a value of 1 being completely similar. Its syntax is:

```julia
compare(::AbstractString, ::AbstractString, ::PreMetric = TokenMax(Levenshtein()))
```

## Distances
- Edit Distances
	- [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) `Hamming()`
	- [Jaro Distance](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) `Jaro()`
	- [Levenshtein Distance](https://en.wikipedia.org/wiki/Levenshtein_distance) `Levenshtein()`
	- [Damerau-Levenshtein Distance](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance) `DamerauLevenshtein()`
	- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html) `RatcliffObershelp()`
- Q-gram distances compare the set of all substrings of length `q` in each string.
	- QGram Distance `Qgram(q)`
	- [Cosine Distance](https://en.wikipedia.org/wiki/Cosine_similarity) `Cosine(q)`
	- [Jaccard Distance](https://en.wikipedia.org/wiki/Jaccard_index) `Jaccard(q)`
	- [Overlap Distance](https://en.wikipedia.org/wiki/Overlap_coefficient) `Overlap(q)`
	- [Sorensen-Dice Distance](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient) `SorensenDice(q)`

- The package includes distance "modifiers", that can be applied to any distance.

	- [Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) boosts the similary score of strings with common prefixes.  The Winkler adjustment was originally defined for the Jaro similarity score but this package defines it for any string distance.
	- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) returns the maximal similarity score between the shorter string and substrings of the longer string.
	- [TokenSort](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders by reording words alphabetically. 
	- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers by comparing the intersection of two strings with each string.
	- [TokenMax](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) combines scores using the base distance, the `Partial`, `TokenSort` and `TokenSet` modifiers, with penalty terms depending on string lengths.

	```julia
	compare("martha", "marhta", Jaro())
	compare("martha", "marhta", Winkler(Jaro()))
	compare("william", "williams", QGram(2))
	compare("william", "williams", Winkler(QGram(2)))
	compare("New York Yankees", "Yankees", Levenshtein())
	compare("New York Yankees", "Yankees", Partial(Levenshtein()))
	compare("mariners vs angels", "los angeles angels at seattle mariners", Jaro())
	compare("mariners vs angels", "los angeles angels at seattle mariners", TokenSet(Jaro()))
	compare("mariners vs angels", "los angeles angels at seattle mariners", TokenMax(RatcliffObershelp()))
	```

## Find
`find_best` returns the index of the element with the highest similarity score.
It returns nothing if all elements have a similarity score below `min_score` (default to 0.0)
```julia
find_best("New York", ["NewYork", "Newark", "San Francisco"], Levenshtein())
#> 1
```

`find_all` returns the indices of the elements with a similarity score higher than a minimum value (default to 0.8)

```julia
find_all("New York", ["NewYork", "Newark", "San Francisco"], Levenshtein(); min_score = 0.8)
#> 1-element Array{String,1}:
#> [1]
```

While these functions are defined for any distance, they are particularly optimized for `Levenshtein` and `DamerauLevenshtein` distances (as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`)

## Evaluate

The function `compare` returns a similarity score: a value of 0 means completely different and a value of 1 means completely similar. In contrast, the function `evaluate` returns the litteral distance between two strings, with a value of 0 being completely similar. Some distances are between 0 and 1, while others are unbouded.

```julia
compare("New York", "New York", Levenshtein())
#> 1.0
evaluate(Levenshtein(), "New York", "New York")
#> 0
```

## Which distance should I use?

As a rule of thumb, 
- Standardize strings before comparing them (cases, whitespaces, accents, abbreviations...)
- The distance `TokenMax(Levenshtein())` (the default for `compare`) is a good choice to link sequence of words (adresses, names) across datasets (see [`fuzzywuzzy`](https://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/))

## References
- [The stringdist Package for Approximate String Matching](https://journal.r-project.org/archive/2014-1/loo.pdf) Mark P.J. van der Loo
- [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)


