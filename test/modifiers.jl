
using StringDistances, Unicode, Random, Test

@testset "Modifiers" begin
	# Partial
	@test Partial(QGram(2))("martha", "marhta") == 6
	@test Partial(QGram(2))("martha", missing) === missing
	@test Partial(Levenshtein())("martha", "marhta") == 2
	@test Partial(RatcliffObershelp())("martha", "marhta") ≈ 0.16666666 atol = 1e-5
	@test Partial(RatcliffObershelp())("martha", "marhtaXXX") ≈ 0.16666666 atol = 1e-5
	@test Partial(RatcliffObershelp())("martha", missing) === missing

	# TokenSort
	@test TokenSort(QGram(2))("martha", "marhta") == 6
	@test TokenSort(QGram(2))("martha", missing) === missing
	@test TokenSort(Levenshtein())("martha", "marhta") == 2
	@test TokenSort(RatcliffObershelp())("martha", "marhta") ≈ 0.16666666 atol = 1e-5

	# TokenSet
	@test TokenSet(QGram(2))("martha", "marhta") == 6
	@test TokenSet(QGram(2))("martha", missing) === missing
	@test TokenSet(Levenshtein())("martha", "marhta") == 2
	@test TokenSet(RatcliffObershelp())("martha", "marhta") ≈ 0.16666666 atol = 1e-5

	# TokenMax
	@test TokenMax(QGram(2))("martha", "marhta") ≈ 0.6
	@test TokenMax(QGram(2))("martha", missing) === missing
	@test TokenMax(Levenshtein())("martha", "marhta") ≈ 1/3
	@test TokenMax(RatcliffObershelp())("martha", "marhta") ≈ 0.16666666 atol = 1e-5
end

@testset "Similarity" begin
	# Qgram
	@test similarity("", "abc", QGram(1)) ≈ 0.0 atol = 1e-4
	@test similarity("abc", "cba", QGram(1)) ≈ 1.0 atol = 1e-4
	@test similarity("abc", "ccc", QGram(1)) ≈ 1/3 atol = 1e-4
	similarity("aüa", "aua", TokenMax(QGram(2)))
	@test similarity("", "abc", Jaccard(2)) ≈ 0.0 atol = 1e-4
	@test similarity("martha", "martha", Jaccard(2)) ≈ 1.0 atol = 1e-4
	@test similarity("martha", "martha", Jaccard(2)) ≈ 1.0 atol = 1e-4
	@test similarity("aa", "aa ", Partial(Jaccard(2))) ≈ 1.0
	@test similarity("martha", "martha", Cosine(2)) ≈ 1.0 atol = 1e-4
	@test similarity("martha", "martha", Overlap(2)) ≈ 1.0 atol = 1e-4
	@test similarity("martha", "martha", SorensenDice(2)) ≈ 1.0 atol = 1e-4

	# Jaro
	@test similarity("aüa", "aua", Hamming()) ≈ 2/3
	@test similarity("aüa", "aua", Jaro()) ≈ 0.77777777 atol = 1e-4
	@test similarity("New York Yankees",  "", Partial(Jaro())) ≈ 0.0

	# JaroWinkler
	@test similarity("martha", "marhta", JaroWinkler()) ≈ 0.9611 atol = 1e-4
	@test similarity("dwayne", "duane", JaroWinkler()) ≈ 0.84 atol = 1e-4
	@test similarity("dixon", "dicksonx", JaroWinkler()) ≈ 0.81333 atol = 1e-4
	@test similarity("william", "williams", JaroWinkler()) ≈ 0.975 atol = 1e-4
	@test similarity("", "foo", JaroWinkler()) ≈ 0.0 atol = 1e-4
	@test similarity("a", "a", JaroWinkler()) ≈ 1.0 atol = 1e-4
	@test similarity("abc", "xyz", JaroWinkler()) ≈ 0.0 atol = 1e-4

	#Levenshtein
	similarity("aüa", "aua", Levenshtein())
	@test similarity("ok", missing, Levenshtein()) === missing
	similarity("aüa", "aua", OptimalStringAlignment())
	@test StringDistances.Normalized(Partial(OptimalStringAlignment()))("ab", "cde") == 1.0
	@test similarity("ab", "de", Partial(OptimalStringAlignment())) == 0

	# RatcliffObershelp
	@test similarity("New York Mets vs Atlanta Braves", "", RatcliffObershelp())  ≈ 0.0
	@test round(Int, 100 * similarity("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", RatcliffObershelp())) == 5
	similarity("aüa", "aua", TokenMax(RatcliffObershelp()))
	@test similarity("New York Yankees",  "Yankees", Partial(RatcliffObershelp())) ≈ 1.0
	@test similarity("New York Yankees",  "", Partial(RatcliffObershelp())) ≈ 0.0
	#@test similarity("mariners vs angels", "los angeles angels at seattle mariners", Partial(RatcliffObershelp())) ≈ 0.444444444444
	@test similarity("HSINCHUANG", "SINJHUAN", Partial(RatcliffObershelp())) ≈ 0.875
	@test similarity("HSINCHUANG", "LSINJHUANG DISTRIC", Partial(RatcliffObershelp())) ≈ 0.8
	@test similarity("HSINCHUANG", "SINJHUANG DISTRICT", Partial(RatcliffObershelp())) ≈ 0.8
	@test similarity("HSINCHUANG",  "SINJHUANG", Partial(RatcliffObershelp())) ≈ 0.8888888888888
	@test similarity("New York Mets vs Atlanta Braves", "Atlanta Braves vs New York Mets", TokenSort(RatcliffObershelp())) ≈ 1.0
	@test similarity(graphemes("New York Mets vs Atlanta Braves"), graphemes("Atlanta Braves vs New York Mets"), Partial(RatcliffObershelp()))  ≈ similarity("New York Mets vs Atlanta Braves", "Atlanta Braves vs New York Mets", Partial(RatcliffObershelp()))
	@test similarity("mariners vs angels", "los angeles angels of anaheim at seattle mariners", TokenSet(RatcliffObershelp())) ≈ 1.0 - 0.09090909090909094
	@test similarity("New York Mets vs Atlanta Braves", "", TokenSort(RatcliffObershelp()))  ≈ 0.0
	@test similarity("mariners vs angels", "", TokenSet(RatcliffObershelp())) ≈ 0.0
	@test similarity("mariners vs angels", "los angeles angels at seattle mariners", TokenSet(Partial(RatcliffObershelp()))) ≈ 1.0
	@test similarity("mariners", "mariner", TokenMax(RatcliffObershelp())) ≈ 0.933333333333333
	@test similarity("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", RatcliffObershelp()) ≈ 5 / 100 atol = 1e-2
	@test similarity("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", Partial(RatcliffObershelp())) ≈ 7 / 100 atol = 1e-2
	@test similarity("mariners", "mariner are playing tomorrow", TokenMax(RatcliffObershelp())) ≈ 79 / 100 atol = 1e-2
	@test similarity("mariners", "mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow", Partial(RatcliffObershelp())) ≈ 88 / 100 atol = 1e-2
	@test similarity("mariners", "mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow", TokenSort(RatcliffObershelp())) ≈ 11 / 100 atol = 1e-2
	@test similarity("mariners", "are mariner playing tomorrow", RatcliffObershelp()) ≈ 39 / 100 atol = 1e-2
	@test similarity("mariners", "are mariner playing tomorrow", Partial(RatcliffObershelp())) ≈ 88 / 100 atol = 1e-2
	@test similarity("mariners", "mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow", TokenSet(RatcliffObershelp())) ≈ 39 / 100 atol = 1e-2
	@test similarity("mariners", "mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow", TokenSet(Partial(RatcliffObershelp()))) ≈ 88 / 100 atol = 1e-2
	# not exactly the same because tokenmax has uses the max of rounded tokenset etc
	@test similarity("mariners", "mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow mariner are playing tomorrow", TokenMax(RatcliffObershelp())) ≈ 78.75 / 100 atol = 1e-2

	@test similarity("martha", "marhta", Levenshtein()) ≈ 2 / 3
	@test_deprecated compare("martha", "marhta", Levenshtein())


	# check min
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
	for dist in (Levenshtein, OptimalStringAlignment)
		for i in eachindex(strings)
			if similarity(strings[i]..., dist()) <  1 / 3
				@test similarity(strings[i]..., dist() ; min_score = 1/ 3) ≈ 0.0
			else
				@test similarity(strings[i]..., dist() ; min_score = 1/ 3) ≈ similarity(strings[i]..., dist())
			end
		end
	end
end


@testset "Find*" begin
	# findnearest
	@test findnearest("New York", ["NewYork", "Newark", "San Francisco"], Levenshtein()) == ("NewYork", 1)
	@test findnearest("New York", ["San Francisco", "NewYork", "Newark"], Levenshtein()) == ("NewYork", 2)
	@test findnearest("New York", ["Newark", "San Francisco", "NewYork"], Levenshtein()) == ("NewYork", 3)
	@test findnearest("New York", ["NewYork", missing, "Newark"], Levenshtein()) == ("NewYork", 1)
	@test findnearest("New York", [missing, missing], Levenshtein()) == (nothing, nothing)
	@test findnearest("New York", ["NewYork", "Newark", "San Francisco"], Jaro()) == ("NewYork", 1)
	@test findnearest("New York", ["NewYork", "Newark", "San Francisco"], QGram(2)) == ("NewYork", 1)
	@test findnearest("New York", ["Newark", "San Francisco", "NewYork"], QGram(2)) == ("NewYork", 3)
	@test findnearest("New York", [missing, "NewYork"], QGram(2)) == ("NewYork", 2)

	# findall
	@test findall("New York", ["NewYork", "Newark", "San Francisco"], Levenshtein()) == [1]
	@test findall("New York", ["NewYork", missing, "Newark"], Levenshtein()) == [1]
	@test findall("x", [missing], Levenshtein(); min_score = 0.0) == Int[]
	@test findall("New York", ["NewYork", "Newark", "San Francisco"], Jaro()) == [1, 2]
	@test findall("New York", ["NewYork", "Newark", "San Francisco"], Jaro(); min_score = 0.99) == Int[]
	@test findall("New York", ["NewYork", "Newark", "San Francisco"], QGram(2); min_score = 0.99) == Int[]
	@test findall("New York", [missing, "NewYork"], QGram(2); min_score = 0.7) == [2]



	if VERSION >= v"1.2.0"
		@test findnearest("New York", skipmissing(["NewYork", "Newark", missing]), Levenshtein()) == ("NewYork", 1)
		@test findnearest("New York", skipmissing(Union{AbstractString, Missing}[missing, missing]), Levenshtein()) == (nothing, nothing)
		@test findall("New York", skipmissing(["NewYork", "Newark", missing]), Levenshtein()) == [1]
		@test findall("New York", skipmissing(Union{AbstractString, Missing}[missing, missing]), Levenshtein()) == []
	end


	Random.seed!(2)
	y = map(Random.randstring, rand(5:25,1_000))
	x = Random.randstring(10)
	for dist in (Levenshtein(), OptimalStringAlignment(), QGram(2), Partial(OptimalStringAlignment()), TokenMax(OptimalStringAlignment()))
		result = [similarity(x, y, dist) for y in y]
		@test findnearest(x, y, dist)[2] == findmax(result)[2]
		@test findall(x, y, dist; min_score = 0.4) == findall(result .>= 0.4)
	end




end
