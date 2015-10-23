[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)


# StringDistances

Edit Distances

- [x] Hamming Distance
- [x] Jaro Distance
- [x] Jaro-Winkler Distance
- [x] Levenshtein Distance
- [x] Damerau-Levenshtein Distance

Q-gram Distances

- [x] qgram
- [x] cosine
- [x] jaccard

Type supports

- [x] ASCIIString
- [x] UTF8String
- [ ] Unicode



Examples
```julia
using StringDistances
hamming("MARTHA", "MARHTA")
levenshtein("MARTHA", "MARHTA")
damerau_levenshtein("MARTHA", "MARHTA")
jaro("MARTHA", "MARHTA")
jaro_winkler("MARTHA", "MARHTA"; scaling_factor = 0.1, boosting_threshold = 0.7, long_threshold = 5)
```

