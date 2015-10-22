
using Base.Test


@test_approx_eq_eps jaro_winkler("MARTHA", "MARHTA", boosting_threshold = 0.0, long_threshold = 100) 0.9611 1e-4
@test_approx_eq_eps jaro_winkler("DWAYNE", "DUANE", boosting_threshold = 0.0, long_threshold = 100) 0.84 1e-4
@test_approx_eq_eps jaro_winkler("DIXON", "DICKSONX", boosting_threshold = 0.0, long_threshold = 100) 0.81333 1e-4
@test_approx_eq_eps jaro_winkler("William", "Williams", boosting_threshold = 0.0, long_threshold = 100) 0.975 1e-4
@test_approx_eq_eps jaro_winkler("", "foo", boosting_threshold = 0.0, long_threshold = 100) 0.0 1e-4
@test_approx_eq_eps jaro_winkler("a", "a", boosting_threshold = 0.0, long_threshold = 100) 1.0 1e-4
@test_approx_eq_eps jaro_winkler("abc", "xyz", boosting_threshold = 0.0, long_threshold = 100) 0.0 1e-4



@test levenshtein("", "") == 0
@test levenshtein("abc", "") == 3
@test levenshtein("", "abc") == 3
@test levenshtein("bc", "abc") == 1
@test levenshtein("kitten", "sitting") == 3
@test levenshtein("Saturday", "Sunday") == 3



@test hamming("", "") == 0
@test hamming("", "abc") == 3
@test hamming("abc", "abc") == 0
@test hamming("acc", "abc") == 1
@test hamming("abcd", "abc") == 1
@test hamming("abc", "abcd") == 1
@test hamming("testing", "this is a test") == 13
@test hamming("Saturday", "Sunday") == 7



