
using StringDistances, Unicode, Test

@testset "Distances" begin

	@testset "Jaro" begin
		@test evaluate(Jaro(), "martha", "marhta") ≈  0.05555555555555547
		@test evaluate(Jaro(), "es an ", " vs an") ≈ 0.2777777777777777
		@test evaluate(Jaro(), " vs an", "es an ") ≈ 0.2777777777777777
		@test evaluate(Jaro(), [1, 2, 3], [1,2, 4]) ≈ 0.2222222222222222
		@test evaluate(Jaro(), graphemes("alborgów"), graphemes("amoniak")) == evaluate(Jaro(), "alborgów", "amoniak")
		@test Jaro()(" vs an", "es an ") ≈ 0.2777777777777777
		@test result_type(Jaro(), "hello", "world") == typeof(float(1))
		@inferred evaluate(Jaro(), "", "")
		@test ismissing(evaluate(Jaro(), "", missing))
	end




	@testset "Levenshtein" begin
		@test evaluate(Levenshtein(), "", "") == 0
		@test evaluate(Levenshtein(), "abc", "") == 3
		@test evaluate(Levenshtein(), "", "abc") == 3
		@test evaluate(Levenshtein(), "bc", "abc") == 1
		@test evaluate(Levenshtein(), "kitten", "sitting") == 3
		@test evaluate(Levenshtein(), "saturday", "sunday") == 3
		@test evaluate(Levenshtein(), "hi, my name is", "my name is") == 4
		@test evaluate(Levenshtein(), "alborgów", "amoniak") == 6
		@test evaluate(Levenshtein(), [1, 2, 3], [1, 2, 4]) == 1
		@test evaluate(Levenshtein(), graphemes("alborgów"), graphemes("amoniak")) == evaluate(Levenshtein(), "alborgów", "amoniak")
		@test Levenshtein()("", "abc") == 3
		@test result_type(Levenshtein(), "hello", "world") == Int
		@inferred evaluate(Levenshtein(), "", "")
		@test ismissing(evaluate(Levenshtein(), "", missing))
	end

	@testset "DamerauLevenshtein" begin
		@test evaluate(DamerauLevenshtein(), "", "") == 0
		@test evaluate(DamerauLevenshtein(), "abc", "") == 3
		@test evaluate(DamerauLevenshtein(), "bc", "abc") == 1
		@test evaluate(DamerauLevenshtein(), "fuor", "four") == 1
		@test evaluate(DamerauLevenshtein(), "abcd", "acb") == 2
		@test evaluate(DamerauLevenshtein(), "cape sand recycling ", "edith ann graham") == 17
		@test evaluate(DamerauLevenshtein(), "jellyifhs", "jellyfish") == 2
		@test evaluate(DamerauLevenshtein(), "ifhs", "fish") == 2
		@test evaluate(DamerauLevenshtein(), [1, 2, 3], [1,2, 4]) == 1
		@test evaluate(DamerauLevenshtein(), graphemes("alborgów"), graphemes("amoniak")) == evaluate(DamerauLevenshtein(), "alborgów", "amoniak")
		@test DamerauLevenshtein()("bc", "abc") == 1
		@test result_type(DamerauLevenshtein(), "hello", "world") == Int
		@inferred evaluate(DamerauLevenshtein(), "", "")
		@test ismissing(evaluate(DamerauLevenshtein(), "", missing))
	end

	@testset "RatcliffObershelp" begin
		@test evaluate(RatcliffObershelp(), "dixon", "dicksonx") ≈ 1 - 0.6153846153846154
		@test evaluate(RatcliffObershelp(), "alexandre", "aleksander") ≈ 1 - 0.7368421052631579
		@test evaluate(RatcliffObershelp(), "pennsylvania",  "pencilvaneya") ≈ 1 - 0.6666666666666
		@test evaluate(RatcliffObershelp(), "",  "pencilvaneya") ≈ 1.0
		@test evaluate(RatcliffObershelp(),"NEW YORK METS", "NEW YORK MEATS") ≈ 1 -  0.962962962963
		@test evaluate(RatcliffObershelp(), "Yankees",  "New York Yankees") ≈ 0.3913043478260869
		@test evaluate(RatcliffObershelp(), "New York Mets",  "New York Yankees") ≈ 0.24137931034482762
		@test evaluate(RatcliffObershelp(), [1, 2, 3], [1,2, 4]) ≈ 1/3
		@test evaluate(RatcliffObershelp(), graphemes("alborgów"), graphemes("amoniak")) == evaluate(RatcliffObershelp(), "alborgów", "amoniak")
		@test RatcliffObershelp()("pennsylvania",  "pencilvaneya") ≈ 1 - 0.6666666666666
		@test result_type(RatcliffObershelp(), "hello", "world") == typeof(float(1))
		@inferred evaluate(RatcliffObershelp(), "", "")
		@test ismissing(evaluate(RatcliffObershelp(), "", missing))
	end


	@testset "QGram" begin
		@test evaluate(QGram(1), "abc", "abc") == 0
		@test evaluate(QGram(1), "", "abc") == 3
		@test evaluate(QGram(1), "abc", "cba") == 0
		@test evaluate(QGram(1), "abc", "ccc") == 4
		@test evaluate(QGram(4), "aü☃", "aüaüafs") == 4
		@test evaluate(QGram(2), SubString("aü☃", 1, 4), SubString("aüaüafs", 1, 4)) == 2
		@test evaluate(QGram(2), graphemes("alborgów"), graphemes("amoniak")) ≈ evaluate(QGram(2), "alborgów", "amoniak")
		@test QGram(1)("abc", "cba") == 0
		@test result_type(QGram(1), "hello", "world") == Int
		@test ismissing(evaluate(QGram(1), "", missing))
		@inferred evaluate(QGram(1), "", "")
	end



	@testset "Cosine" begin
		@test isnan(evaluate(Cosine(2), "", "abc"))
		@test evaluate(Cosine(2), "abc", "ccc") ≈ 1 atol = 1e-4
		@test evaluate(Cosine(2), "leia", "leela") ≈ 0.7113249 atol = 1e-4
		@test evaluate(Cosine(2), [1, 2, 3], [1, 2, 4]) ≈ 0.5
		@test evaluate(Cosine(2), graphemes("alborgów"), graphemes("amoniak")) ≈ evaluate(Cosine(2), "alborgów", "amoniak")
		@test Cosine(2)("leia", "leela") ≈ 0.7113249 atol = 1e-4
		@test result_type(Cosine(2), "hello", "world") == typeof(float(1))
		@inferred evaluate(Cosine(2), "", "")
		@test ismissing(evaluate(Cosine(2), "", missing))
	end

	@testset "Jaccard" begin
		@test evaluate(Jaccard(1), "", "abc") ≈ 1.0
		@test evaluate(Jaccard(1), "abc", "ccc") ≈ 2/3 atol = 1e-4
		@test evaluate(Jaccard(2), "leia", "leela") ≈ 0.83333 atol = 1e-4
		@test evaluate(Jaccard(2), [1, 2, 3], [1, 2, 4]) ≈ 2/3 atol = 1e-4
		@test evaluate(Jaccard(2), graphemes("alborgów"), graphemes("amoniak")) ≈ evaluate(Jaccard(2), "alborgów", "amoniak")
		@test Jaccard(2)("leia", "leela") ≈ 0.83333 atol = 1e-4
		@test result_type(Jaccard(1), "hello", "world") == typeof(float(1))
		@inferred evaluate(Jaccard(1), "", "")
		@test ismissing(evaluate(Jaccard(1), "", missing))
	end

	@testset "SorensenDice" begin
		@test evaluate(SorensenDice(1), "night", "nacht") ≈ 0.4 atol = 1e-4
		@test evaluate(SorensenDice(2), "night", "nacht") ≈ 0.75 atol = 1e-4
		@test evaluate(SorensenDice(2), graphemes("alborgów"), graphemes("amoniak")) ≈ evaluate(SorensenDice(2), "alborgów", "amoniak")
		@test SorensenDice(2)("night", "nacht") ≈ 0.75 atol = 1e-4
		@test result_type(SorensenDice(1), "hello", "world") == typeof(float(1))
		@inferred evaluate(SorensenDice(1), "", "")
		@test ismissing(evaluate(SorensenDice(1), "", missing))
	end

	@testset "Overlap" begin
		@test evaluate(Overlap(1), "night", "nacht") ≈ 0.4 atol = 1e-4
		@test evaluate(Overlap(1), "context", "contact") ≈ .2 atol = 1e-4
		@test Overlap(1)("context", "contact") ≈ .2 atol = 1e-4
		@test result_type(Overlap(1), "hello", "world") == typeof(float(1))
		@inferred evaluate(Overlap(1), "", "")
		@test ismissing(evaluate(Overlap(1), "", missing))
	end




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
			(Jaro(), [0.05555556 0.17777778 0.23333333 0.04166667 1.00000000 0.00000000 1.00000000 0.44444444 0.25396825 0.2805556 0.2285714 0.48809524 0.3916667 0.07407407 0.16666667 0.21666667]),
			(QGram(1), [0   3   3   1 3  0   6   4   5   4   4  11  14   0   0   3]),
			(QGram(2), [  6   7   7   1 2 0   4   4   7   8   4  13  32   8   6   5]),
			(Jaccard(1), [0.0 0.4285714 0.3750000 0.1666667       1.0 0.0 1.0000000 0.6666667 0.5714286 0.3750000 0.2000000 0.8333333 0.5000000 0.0 0.0 0.2500000]),
			(Jaccard(2),  [ 0.7500000 0.8750000 0.7777778 0.1428571       1.0     NaN 1.0000000 1.0000000 0.7777778 0.8000000 0.3076923 1.0000000 0.9696970 0.6666667 1.0000000 0.8333333]),
			(Cosine(2), [0.6000000 0.7763932 0.6220355 0.0741799  NaN  NaN 1.0000000 1.0000000 0.6348516 0.6619383 0.1679497 1.0000000 0.9407651 0.5000000 1.0000000 0.7113249]))
	# Test with R package StringDist
	for x in solutions
		t, solution = x
		for i in 1:length(solution)
			if isnan(evaluate(t, strings[i]...))
				@test isnan(solution[i])
			else
				@test evaluate(t, strings[i]...) ≈ solution[i] atol = 1e-4
			end
		end
	end
	# test  RatcliffObershelp
	solution = [83, 73, 62, 93, 0, 100, 0, 33, 62, 71, 83, 27, 33, 78, 50, 67]
	for i in eachindex(strings)
		@test round(Int, (1 - evaluate(RatcliffObershelp(), strings[i]...)) * 100) ≈ solution[i] atol = 1e-4
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




#= Fuzzywuzzy usesRatcliffObershelp  if python-Levenshtein not installed, fuzzywuzzy uses RatcliffObershelp)
from fuzzywuzzy import fuzz
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
for x in strings:
   print(fuzz.ratio(x[0], x[1]))
=#

