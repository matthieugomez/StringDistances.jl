
using StringDistances, Test


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