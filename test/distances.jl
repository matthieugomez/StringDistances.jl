
using StringDistances, Base.Test


@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "martha", "marhta") 1 - 0.9611 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "dwayne", "duane") 1 - 0.84 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "dixon", "dicksonx") 1 - 0.81333 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "william", "williams") 1 - 0.975 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "", "foo") 1.0 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "a", "a") 0.0 1e-4
@test_approx_eq_eps evaluate(JaroWinkler(0.1, 0.0, 100), "abc", "xyz") 1.0 1e-4


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


@test_approx_eq_eps evaluate(Jaccard(1), "", "abc") 1.0 1e-4
@test_approx_eq_eps evaluate(Jaccard(1), "abc", "ccc") .666666 1e-4
@test_approx_eq_eps evaluate(Jaccard(2), "leia", "leela") 0.83333 1e-4


@test_approx_eq_eps evaluate(Cosine(2), "", "abc") 1 1e-4
@test_approx_eq_eps evaluate(Cosine(2), "abc", "ccc") 1 1e-4
@test_approx_eq_eps evaluate(Cosine(2), "leia", "leela") 0.7113249 1e-4
