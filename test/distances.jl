
using StringDistances, Base.Test


@test_approx_eq_eps jaro_winkler("MARTHA", "MARHTA", boosting_threshold = 0.0, long_threshold = 100) 1 - 0.9611 1e-4
@test_approx_eq_eps jaro_winkler("DWAYNE", "DUANE", boosting_threshold = 0.0, long_threshold = 100) 1 - 0.84 1e-4
@test_approx_eq_eps jaro_winkler("DIXON", "DICKSONX", boosting_threshold = 0.0, long_threshold = 100) 1 - 0.81333 1e-4
@test_approx_eq_eps jaro_winkler("William", "Williams", boosting_threshold = 0.0, long_threshold = 100) 1 - 0.975 1e-4
@test_approx_eq_eps jaro_winkler("", "foo", boosting_threshold = 0.0, long_threshold = 100) 1.0 1e-4
@test_approx_eq_eps jaro_winkler("a", "a", boosting_threshold = 0.0, long_threshold = 100) 0.0 1e-4
@test_approx_eq_eps jaro_winkler("abc", "xyz", boosting_threshold = 0.0, long_threshold = 100) 1.0 1e-4



@test levenshtein("", "") == 0
@test levenshtein("abc", "") == 3
@test levenshtein("", "abc") == 3
@test levenshtein("bc", "abc") == 1
@test levenshtein("kitten", "sitting") == 3
@test levenshtein("Saturday", "Sunday") == 3


@test damerau_levenshtein("", "") == 0
@test damerau_levenshtein("abc", "") == 3
@test damerau_levenshtein("bc", "abc") == 1
@test damerau_levenshtein("fuor", "four") == 1
@test damerau_levenshtein("abcd", "acb") == 2
@test damerau_levenshtein("cape sand recycling ", "edith ann graham") == 17
@test damerau_levenshtein("jellyifhs", "jellyfish") == 2
@test damerau_levenshtein("ifhs", "fish") == 2


@test hamming("", "") == 0
@test hamming("", "abc") == 3
@test hamming("abc", "abc") == 0
@test hamming("acc", "abc") == 1
@test hamming("abcd", "abc") == 1
@test hamming("abc", "abcd") == 1
@test hamming("testing", "this is a test") == 13
@test hamming("Saturday", "Sunday") == 7


@test qgram("", "abc", q = 1) == 3
@test qgram("abc", "cba", q = 1) == 0
@test qgram("abc", "ccc", q = 1) == 4


@test_approx_eq_eps jaccard("", "abc", q = 1) 1.0 1e-4
@test_approx_eq_eps jaccard("abc", "ccc", q = 1) .666666 1e-4
@test_approx_eq_eps jaccard("leia", "leela", q = 2) 0.83333 1e-4


@test_approx_eq_eps cosine("", "abc", q = 2) 1 1e-4
@test_approx_eq_eps cosine("abc", "ccc", q = 2) 1 1e-4
@test_approx_eq_eps cosine("leia", "leela", q = 2) 0.7113249 1e-4
