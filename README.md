[![Build status](https://github.com/matthieugomez/StringDistances.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/StringDistances.jl/actions)


## Installation
The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and so can be installed at the REPL with `] add StringDistances`.

## Supported Distances
String distances act over any pair of iterators that define `length` (e.g. `AbstractStrings`, `GraphemeIterators`, or `AbstractVectors`)

The available distances are:
- Edit Distances
	- Hamming Distance `Hamming() <: SemiMetric`
	- [Jaro and Jaro-Winkler Distance](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) `Jaro()` `JaroWinkler() <: SemiMetric`
	- [Levenshtein Distance](https://en.wikipedia.org/wiki/Levenshtein_distance) `Levenshtein() <: Metric`
	- [Optimal String Alignment Distance](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance#Optimal_string_alignment_distance) (a.k.a. restricted Damerau-Levenshtein) `OptimalStringAlignment() <: SemiMetric`
	- [Damerau-Levenshtein Distance](https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance#Distance_with_adjacent_transpositions) `DamerauLevenshtein() <: Metric`
	- [RatcliffObershelp Distance](https://xlinux.nist.gov/dads/HTML/ratcliffObershelp.html) `RatcliffObershelp() <: SemiMetric`
- Q-gram distances compare the set of all substrings of length `q` in each string (and which 
	- QGram Distance `Qgram(q::Int) <: SemiMetric`
	- [Cosine Distance](https://en.wikipedia.org/wiki/Cosine_similarity) `Cosine(q::Int) <: SemiMetric`
	- [Jaccard Distance](https://en.wikipedia.org/wiki/Jaccard_index) `Jaccard(q::Int) <: SemiMetric`
	- [Overlap Distance](https://en.wikipedia.org/wiki/Overlap_coefficient) `Overlap(q::Int) <: SemiMetric`
	- [Sorensen-Dice Distance](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient) `SorensenDice(q::Int) <: SemiMetric`
	- [MorisitaOverlap Distance](https://en.wikipedia.org/wiki/Morisita%27s_overlap_index) `MorisitaOverlap(q::Int) <: SemiMetric`
	- [Normalized Multiset Distance](https://www.sciencedirect.com/science/article/pii/S1047320313001417) `NMD(q::Int) <: SemiMetric`

## Syntax
Following the `Distances.jl` package, string distances can inherit from two abstract types: `StringSemiMetric <: SemiMetric` or `StringMetric <: Metric`.
## Computing the distance between two strings (or iterators)
You can always compute a certain distance between two strings  using the following syntax
```julia
r = evaluate(dist, x, y)
r = dist(x, y)
```
Here, `dist` is an instance of a distance type: for example, the type for the Levenshtein distance is `Levenshtein`. You can compute the Levenshtein distance between `x` and `y` as
```julia
r = evaluate(Levenshtein(), x, y)
r = Levenshtein()(x, y)
```

The function `compare` returns the similarity score, defined as 1 minus the normalized distance between two strings. It always returns an element of type `Float64`. A value of 0.0 means completely different and a value of 1.0 means completely similar.

```julia
Levenshtein()("martha", "martha")
#> 0
compare("martha", "martha", Levenshtein())
#> 1.0
```

## Computing the distance between two AbstractVectors of strings (or iterators)
Consider `X` and `Y` two `AbstractVectors` of iterators. You can compute the matrix of distances across elements, `dist(X[i], Y[j])`, as follows:
```julia
pairwise(dist, X, Y)
```

For instance, 
```julia
pairwise(Jaccard(3), ["martha", "kitten"], ["marhta", "sitting"])
```
`pairwise` is optimized in various ways (e.g., for the case of QGram-distances, dictionary of qgrams are pre-computed)

## Find closest string
The package also adds convenience functions to find elements in a iterator of strings closest to a given string

- `findnearest` returns the value and index of the element in `itr` with the highest similarity score with `s`. Its syntax is:
	```julia
	findnearest(s, itr, dist)
	```

- `findall` returns the indices of all elements in `itr` with a similarity score with `s` higher than a minimum score. Its syntax is:
	```julia
	findall(s, itr, dist; min_score = 0.8)
	```

The functions `findnearest` and `findall` are particularly optimized for the `Levenshtein` and `OptimalStringAlignment` distances, as these algorithm can stop early if the distance becomes higher than a certain threshold.



### fuzzywuzzy
The package also defines Distance "modifiers" that are inspired by the Python package - [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/). These modifiers are particularly helpful to match strings composed of multiple words (e.g. addresses, company names).
- [Partial](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) returns the minimum of the distance between the shorter string and substrings of the longer string.
- [TokenSort](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders by returning the distance of the two strings, after re-ordering words alphabetically. 
- [TokenSet](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) adjusts for differences in word orders and word numbers by returning the distance between the intersection of two strings with each string.
- [TokenMax](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/) normalizes the distance, and combine the `Partial`, `TokenSort` and `TokenSet` modifiers, with penalty terms depending on string.   `TokenMax(Levenshtein())` corresponds to the distance defined in [fuzzywuzzy](http://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)


```julia
Levenshtein()("this string", "this string is longer") = 10
Partial(Levenshtein())("this string", "this string is longer") = 0
```



## Notes
- All string distances are case sensitive.


