[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)


# StringDistances

ASCII

- [x] Hamming Distance
- [x] Jaro Distance and Jaro-Winkler Distance
- [x] Levenshtein Distance
- [x] Damerau-Levenshtein Distance
- [x] Qgram Distance
- [x] Cosine Distance
- [x] Jaccard Distance

AbstractString

- [x] Hamming Distance
- [] Jaro Distance and Jaro-Winkler Distance
- [x] Levenshtein Distance
- [x] Damerau-Levenshtein Distance
- [x] Qgram Distance
- [x] Cosine Distance
- [x] Jaccard Distance



Examples
```julia
using StringDistances
hamming("MARTHA", "MARHTA")
levenshtein("MARTHA", "MARHTA")
damerau_levenshtein("MARTHA", "MARHTA")
jaro("MARTHA", "MARHTA")
jaro_winkler("MARTHA", "MARHTA"; scaling_factor = 0.1, boosting_threshold = 0.7, long_threshold = 5)
```

