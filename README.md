[![Build Status](https://travis-ci.org/matthieugomez/StringDistances.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/StringDistances.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/StringDistances.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/StringDistances.jl?branch=master)


# StringDistances

- [x] Hamming Distance
- [x] Jaro Distance and Jaro-Winkler Distance
- [x] Levenshtein Distance
- [x] Damerau-Levenshtein Distance
- [x] Qgram Distance
- [x] Cosine Distance
- [x] Jaccard Distance

Support for ASCII, UTF-8 and Unicode

# Syntax
There are two possible syntaxes for each distance:
```julia
using StringDistances
evaluate(Jaccard(2), "martha", "marhta")
jaccard("martha", "marhta"; q = 2)
```

