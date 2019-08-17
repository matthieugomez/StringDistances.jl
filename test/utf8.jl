using StringDistances, Test

# check with weird utf8 strings
compare("aüa", "aua", TokenMax(RatcliffObershelp()))
compare("aüa", "aua", TokenMax(QGram(2)))
compare("aüa", "aua", DamerauLevenshtein())
compare("aüa", "aua", Hamming())
compare("aüa", "aua", Jaro())
compare("aüa", "aua", Levenshtein())


s1 = "aü☃"
s2 = "aüaüafs"
dist = QGram(4)
@test evaluate(dist, s1, s2) == 4

# check Substrings work
s1 = SubString(s1, 1, 4)
s2 = SubString(s2, 1, 4)
dist = QGram(2)
@test evaluate(dist, s1, s2) == 2
