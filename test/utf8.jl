using StringDistances, Test

# check with weird utf8 strings
compare(TokenMax(RatcliffObershelp()), "aüa", "aua")
compare(TokenMax(QGram(2)), "aüa", "aua")
compare(DamerauLevenshtein(), "aüa", "aua")
compare(Hamming(), "aüa", "aua")
compare(Jaro(), "aüa", "aua")
compare(Levenshtein(), "aüa", "aua")


s1 = "aü☃"
s2 = "aüaüafs"
dist = QGram(4)
@test evaluate(dist, s1, s2) == 4

# check Substrings work
s1 = SubString(s1, 1, 4)
s2 = SubString(s2, 1, 4)
dist = QGram(2)
@test evaluate(dist, s1, s2) == 2
