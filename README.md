[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)

## Installation
The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and so can be installed at the REPL with `] add StringDistances`.

## Evaluate
The function `evaluate` returns the distance between two strings. Its syntax is:

```julia
evaluate(dist, s1, s2)
```

where `s1` and `s2` can be any iterator with a `length` method (e.g. `AbstractString`, `GraphemeIterator`, `AbstractVector`...), and `dist` is one of the following distances:

- Edit Distances
	- [Jaro Distance](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) `Jaro()`
	- [Levenshtein Distance](https://en.wikipedia.org/wiki/Levenshtein_distance) `Levenshtein()`
	- [Damerau-Levenshtein Distance](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance) `DamerauLevenshtein()`
	- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html) `RatcliffObershelp()`
- Q-gram distances compare the set of all substrings of length `q` in each string.
	- QGram Distance `Qgram(q::Int)`
	- [Cosine Distance](https://en.wikipedia.org/wiki/Cosine_similarity) `Cosine(q::Int)`
	- [Jaccard Distance](https://en.wikipedia.org/wiki/Jaccard_index) `Jaccard(q::Int)`
	- [Overlap Distance](https://en.wikipedia.org/wiki/Overlap_coefficient) `Overlap(q::Int)`
	- [Sorensen-Dice Distance](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient) `SorensenDice(q::Int)`

- The package includes distance "modifiers", that can be applied to any distance.

	- [Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) diminishes the distance of strings with common prefixes.  The Winkler adjustment was originally defined for the Jaro similarity score but this package defines it for any string distance.
	- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) returns the minimum distance between the shorter string and substrings of the longer string.
	- [TokenSort](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders by reording words alphabetically. 
	- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers by comparing the intersection of two strings with each string.
	- [TokenMax](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) combines scores using the base distance, the `Partial`, `TokenSort` and `TokenSet` modifiers, with penalty terms depending on string lengths.

Some examples:
```julia
evaluate(Jaro(), "martha", "marhta")
evaluate(Winkler(Jaro()), "martha", "marhta")
evaluate(QGram(2), "martha", "marhta")
evaluate(Winkler(QGram(2)), "martha", "marhta")
evaluate(Levenshtein(), "martha", "marhta")
evaluate(Partial(Levenshtein()), "martha", "marhta")
evaluate(Jaro(), "martha", "marhta")
evaluate(TokenSet(Jaro()), "martha", "marhta")
evaluate(TokenMax(RatcliffObershelp()), "martha", "marhta")
```

Alternatively, each distance can be used as a callable to call the evaluate function of each metric or modified metric, for example:
```julia
Jaro()("martha", "marhta")
Winkler(Jaro())("martha", "marhta")
QGram(2)("martha", "marhta")
```

A good distance to match strings composed of multiple words (like addresses) is `TokenMax(Levenshtein())` (see [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)).

## Compare
The function `compare` is defined as 1 minus the normalized distance between two strings. It always returns a number between 0 and 1: a value of 0 means completely different and a value of 1 means completely similar.
```julia
evaluate(Levenshtein(), "New York", "New York")
#> 0
compare("New York", "New York", Levenshtein())
#> 1.0

```


## Find
- `findmax` returns the value and index of the element in `itr` with the highest similarity score with `s`. Its syntax is:
	```julia
	findmax(s, itr, dist::StringDistance; min_score = 0.0)
	```

- `findall` returns the indices of all elements in `itr` with a similarity score with `s` higher than a minimum value (default to 0.8). Its syntax is:
	```julia
	findall(s, itr, dist::StringDistance; min_score = 0.8)
	```

The functions `findmax` and `findall` are particularly optimized for `Levenshtein` and `DamerauLevenshtein` distances (as well as their modifications via `Partial`, `TokenSort`, `TokenSet`, or `TokenMax`).


## References
- [The stringdist Package for Approximate String Matching](https://journal.r-project.org/archive/2014-1/loo.pdf) Mark P.J. van der Loo
- [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)


