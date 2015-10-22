# StringDistances
String Distances in Julia

- [x] Hamming Distance
- [x] Jaro distance
- [x] Jaro-Winkler Distance
- [x] Levenshtein distance
- [] Damerau-Levenshtein Distance
- [] qgram


Examples
```julia
using StringDistances
hamming("MARTHA", "MARHTA")
levenshtein("MARTHA", "MARHTA")
jaro("MARTHA", "MARHTA")
jaro_winkler("MARTHA", "MARHTA"; scaling_factor = 0.1, boosting_threshold = 0.7, long_threshold = 5)
```

