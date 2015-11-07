
using StringDistances, Base.Test


@test evaluate(Levenshtein(), "", "") == 0
@test evaluate(Levenshtein(), "abc", "") == 3
@test evaluate(Levenshtein(), "", "abc") == 3
@test evaluate(Levenshtein(), "bc", "abc") == 1
@test evaluate(Levenshtein(), "kitten", "sitting") == 3
@test evaluate(Levenshtein(), "saturday", "sunday") == 3

@test evaluate(Levenshtein(), "hi, my name is", "my name is") == 4
@test evaluate(Levenshtein(), "alborgów", "amoniak") == 6

@test evaluate(DamerauLevenshtein(), "", "") == 0
@test evaluate(DamerauLevenshtein(), "abc", "") == 3
@test evaluate(DamerauLevenshtein(), "bc", "abc") == 1
@test evaluate(DamerauLevenshtein(), "fuor", "four") == 1
@test evaluate(DamerauLevenshtein(), "abcd", "acb") == 2
@test evaluate(DamerauLevenshtein(), "cape sand recycling ", "edith ann graham") == 17
@test evaluate(DamerauLevenshtein(), "jellyifhs", "jellyfish") == 2
@test evaluate(DamerauLevenshtein(), "ifhs", "fish") == 2

@test evaluate(Hamming(), "", "") == 0
@test evaluate(Hamming(), "", "abc") == 3
@test evaluate(Hamming(), "abc", "abc") == 0
@test evaluate(Hamming(), "acc", "abc") == 1
@test evaluate(Hamming(), "abcd", "abc") == 1
@test evaluate(Hamming(), "abc", "abcd") == 1
@test evaluate(Hamming(), "testing", "this is a test") == 13
@test evaluate(Hamming(), "saturday", "sunday") == 7

@test evaluate(QGram(1), "", "abc") == 3
@test evaluate(QGram(1), "abc", "cba") == 0
@test evaluate(QGram(1), "abc", "ccc") == 4
@test isnan(evaluate(Cosine(2), "", "abc"))
@test_approx_eq_eps evaluate(Cosine(2), "abc", "ccc") 1 1e-4
@test_approx_eq_eps evaluate(Cosine(2), "leia", "leela") 0.7113249 1e-4
@test_approx_eq evaluate(Jaccard(1), "", "abc") 1.0
@test_approx_eq_eps evaluate(Jaccard(1), "abc", "ccc") .666666 1e-4
@test_approx_eq_eps evaluate(Jaccard(2), "leia", "leela") 0.83333 1e-4
@test_approx_eq_eps evaluate(SorensenDice(1), "night", "nacht") 0.4 1e-4
@test_approx_eq_eps evaluate(SorensenDice(2), "night", "nacht") 0.75 1e-4
@test_approx_eq_eps evaluate(Overlap(1), "night", "nacht") 0.4 1e-4
@test_approx_eq_eps evaluate(Overlap(1), "context", "contact") .2 1e-4

@test matching_blocks("dwayne", "duane") ==
Set([(5,4,2) (1,1,1) (3,3,1)])
@test matching_blocks("dixon", "dicksonx") ==
Set([(1,1,2) (4,6,2)])
@test_approx_eq evaluate(RatcliffObershelp(), "dixon", "dicksonx") 1 - 0.6153846153846154
@test_approx_eq evaluate(RatcliffObershelp(), "alexandre", "aleksander") 1 - 0.7368421052631579
@test_approx_eq evaluate(RatcliffObershelp(), "pennsylvania",  "pencilvaneya") 1 - 0.6666666666666
@test_approx_eq evaluate(RatcliffObershelp(), "",  "pencilvaneya") 1.0
@test_approx_eq evaluate(RatcliffObershelp(),"NEW YORK METS", "NEW YORK MEATS") 1 -  0.962962962963
@test_approx_eq evaluate(RatcliffObershelp(), "Yankees",  "New York Yankees") 0.3913043478260869
@test_approx_eq evaluate(RatcliffObershelp(), "New York Mets",  "New York Yankees") 0.24137931034482762




strings = [
("martha", "marhta"),
("dwayne", "duane") ,
("dixon", "dicksonx"),
("william", "williams"),
("", "foo"),
("a", "a"),
("abc", "xyz"),
("abc", "ccc"),
("kitten", "sitting"),
("saturday", "sunday"),
("hi, my name is", "my name is"),
("alborgów", "amoniak"),
("cape sand recycling ", "edith ann graham"),
( "jellyifhs", "jellyfish"),
("ifhs", "fish"),
("leia", "leela"),
]

solutions = ((Levenshtein(), [2  2  4  1  3  0  3  2  3  3  4  6 17  3  3  2]),
		(DamerauLevenshtein(), [1  2  4  1  3  0  3  2  3  3  4  6 17  2  2  2]),
		(Jaro(), [0.05555556 0.17777778 0.23333333 0.04166667 1.00000000 0.00000000 1.00000000 0.44444444 0.25396825 0.24722222 0.16190476 0.48809524 0.49166667 0.07407407 0.16666667 0.21666667]),
		(QGram(1), [0   3   3   1 3  0   6   4   5   4   4  11  14   0   0   3]),
		(QGram(2), [  6   7   7   1 2 0   4   4   7   8   4  13  32   8   6   5]),
		(Jaccard(1), [0.0 0.4285714 0.3750000 0.1666667       1.0 0.0 1.0000000 0.6666667 0.5714286 0.3750000 0.2000000 0.8333333 0.5000000 0.0 0.0 0.2500000]),
		(Jaccard(2),  [ 0.7500000 0.8750000 0.7777778 0.1428571       1.0     NaN 1.0000000 1.0000000 0.7777778 0.8000000 0.3076923 1.0000000 0.9696970 0.6666667 1.0000000 0.8333333]),
		(Cosine(2), [0.6000000 0.7763932 0.6220355 0.0741799  NaN  NaN 1.0000000 1.0000000 0.6348516 0.6619383 0.1679497 1.0000000 0.9407651 0.5000000 1.0000000 0.7113249]))
# Test with R package StringDist
for x in solutions
	t, solution = x
	for i in 1:length(solution)
		@test_approx_eq_eps evaluate(t, strings[i]...) solution[i] 1e-4
	end
end


#= R test
library(stringdist)
strings = matrix(data = c(
"martha", "marhta",
"dwayne", "duane",
"dixon", "dicksonx",
"william", "williams",
"", "foo",
"a", "a",
"abc", "xyz",
"abc", "ccc",
"kitten", "sitting",
"saturday", "sunday",
"hi, my name is", "my name is",
"alborgów", "amoniak",
"cape sand recycling ", "edith ann graham",
 "jellyifhs", "jellyfish",
"ifhs", "fish",
"leia", "leela"), 
nrow = 2
)
stringdist(strings[1,], strings[2,], method = "jw", p = 0)
stringdist(strings[1,], strings[2,], method = "jw", p = 0.1)
stringdist(strings[1,], strings[2,], method = "qgram", q = 1)

=#


# grapheme
strings = [
(graphemes("martha"), graphemes("marhta")),
(graphemes("dwayne"), graphemes("duane") ),
(graphemes("dixon"), graphemes("dicksonx")),
(graphemes("william"), graphemes("williams")),
(graphemes(""), graphemes("foo")),
(graphemes("a"), graphemes("a")),
(graphemes("abc"), graphemes("xyz")),
(graphemes("abc"), graphemes("ccc")),
(graphemes("kitten"), graphemes("sitting")),
(graphemes("saturday"), graphemes("sunday")),
(graphemes("hi, my name is"), graphemes("my name is")),
(graphemes("alborgów"), graphemes("amoniak")),
(graphemes("cape sand recycling "), graphemes("edith ann graham")),
(graphemes( "jellyifhs"), graphemes("jellyfish")),
(graphemes("ifhs"), graphemes("fish")),
(graphemes("leia"), graphemes("leela")),
]


for x in solutions
	t, solution = x
	for i in 1:length(solution)
		@test_approx_eq_eps evaluate(t, strings[i]...) solution[i] 1e-4
	end
end

@test evaluate(Hamming(), "b\u0300", "a") == 2
@test evaluate(Hamming(), graphemes("b\u0300"), graphemes("a")) == 1