
using StringDistances, Test

# Compare
@test compare("", "abc", Hamming()) ≈ 0.0 atol = 1e-4
@test compare("acc", "abc", Hamming()) ≈ 2/3 atol = 1e-4
@test compare("saturday", "sunday", Hamming()) ≈ 1/8 atol = 1e-4

@test compare("", "abc", QGram(1)) ≈ 0.0 atol = 1e-4
@test compare("abc", "cba", QGram(1)) ≈ 1.0 atol = 1e-4
@test compare("abc", "ccc", QGram(1)) ≈ 1/3 atol = 1e-4

@test compare("", "abc", Jaccard(2)) ≈ 0.0 atol = 1e-4

@test compare("martha", "martha", Jaccard(2)) ≈ 1.0 atol = 1e-4
@test compare("martha", "martha", Cosine(2)) ≈ 1.0 atol = 1e-4
@test compare("martha", "martha", Jaccard(2)) ≈ 1.0 atol = 1e-4
@test compare("martha", "martha", Overlap(2)) ≈ 1.0 atol = 1e-4
@test compare("martha", "martha", SorensenDice(2)) ≈ 1.0 atol = 1e-4


# Winkler
@test compare("martha", "marhta", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.9611 atol = 1e-4
@test compare("dwayne", "duane", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.84 atol = 1e-4
@test compare("dixon", "dicksonx", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.81333 atol = 1e-4
@test compare("william", "williams", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.975 atol = 1e-4
@test compare("", "foo", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.0 atol = 1e-4
@test compare("a", "a", Winkler(Jaro(), 0.1, 0.0)) ≈ 1.0 atol = 1e-4
@test compare("abc", "xyz", Winkler(Jaro(), 0.1, 0.0)) ≈ 0.0 atol = 1e-4

strings = [
("martha", "marhta"),
("dwayne", "duane") ,
("dixon", "dicksonx"),
("william", "williams"),
("", "foo")
]
solutions = [0.03888889 0.16000000 0.18666667 0.02500000 1.00000000]
for i in 1:length(solutions)
	@test compare(strings[i]..., Winkler(Jaro(), 0.1, 0.0)) ≈ (1 - solutions[i]) atol = 1e-4
end





# Partial
@test compare("aa", "aa ", Partial(Jaccard(2))) ≈ 1.0

@test compare("New York Yankees",  "Yankees", Partial(RatcliffObershelp())) ≈ 1.0
@test compare("New York Yankees",  "", Partial(RatcliffObershelp())) ≈ 0.0
@test compare("mariners vs angels", "los angeles angels at seattle mariners", Partial(RatcliffObershelp())) ≈ 0.444444444444


s = "HSINCHUANG"
@test compare(s, "SINJHUAN", Partial(RatcliffObershelp())) ≈ 0.875
@test compare(s, "LSINJHUANG DISTRIC", Partial(RatcliffObershelp())) ≈ 0.8
@test compare(s, "SINJHUANG DISTRICT", Partial(RatcliffObershelp())) ≈ 0.8
@test compare(s,  "SINJHUANG", Partial(RatcliffObershelp())) ≈ 0.8888888888888

@test compare("New York Yankees",  "Yankees", Partial(Hamming())) ≈ 1
@test compare("New York Yankees",  "", Partial(Hamming())) ≈ 1



# Token
@test compare("New York Mets vs Atlanta Braves", "Atlanta Braves vs New York Mets", TokenSort(RatcliffObershelp()))  ≈ 1.0
@test compare("mariners vs angels", "los angeles angels of anaheim at seattle mariners", TokenSet(RatcliffObershelp())) ≈ 1.0 - 0.09090909090909094


@test compare("New York Mets vs Atlanta Braves", "", TokenSort(RatcliffObershelp()))  ≈ 0.0
@test compare("mariners vs angels", "", TokenSet(RatcliffObershelp())) ≈ 0.0




@test compare("mariners", "mariner", TokenMax(RatcliffObershelp())) ≈ 0.933333333333333

#@test_approx_eq compare(TokenSort(RatcliffObershelp()), graphemeiterator("New York Mets vs Atlanta Braves"), graphemeiterator("Atlanta Braves vs New York Mets"))  1.0
#@test_approx_eq compare(TokenSet(RatcliffObershelp()),graphemeiterator("mariners vs angels"), graphemeiterator("los angeles angels of anaheim at seattle mariners")) 1.0 - 0.09090909090909094
#@test_approx_eq compare(TokenSort(RatcliffObershelp()), graphemeiterator("New York Mets vs Atlanta Braves"), graphemeiterator(""))  0.0
#@test_approx_eq compare(TokenSet(RatcliffObershelp()),graphemeiterator("mariners vs angels"), graphemeiterator("")) 0.0



@test compare("mariners vs angels", "los angeles angels at seattle mariners", TokenSet(Partial(RatcliffObershelp()))) ≈ 1.0


@test round(Int, 100 * compare("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", RatcliffObershelp())) == 5


# test with fuzz ratio
@test round(Int, 100 * compare("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", RatcliffObershelp())) == 5
@test round(Int, 100 * compare("为人子女者要堂堂正正做人，千万不可作奸犯科，致使父母蒙羞", "此前稍早些时候中国商务部发布消息称，中美经贸高级别磋商双方牵头人通话，中方就美拟9月1日加征关税进行了严正交涉。", Partial(RatcliffObershelp()))) == 7


