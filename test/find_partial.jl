using StringDistances: get_max_dist, update_max_dist

@testset "`find*_partial`" begin
    @testset "`get_max_dist` and `update_max_dist`" begin
        d = 2
        @testset "$T" for T in (DamerauLevenshtein, Levenshtein, Hamming)
            dist = T(d)
            @test get_max_dist(dist) == d
            @test update_max_dist(dist, 2*d) == T(2*d)

            d1 = 0.5
            @testset "$mod" for mod in (StringDistances.Normalized, StringDistances.TokenMax)
                mod_dist = mod(dist, d1)
                @test get_max_dist(mod_dist) == d1
                @test update_max_dist(mod_dist, 2*d1) == mod(dist, 2*d1)
            end

            @testset "$mod" for mod in (Partial, StringDistances.TokenSort, StringDistances.TokenSet)
                if mod == StringDistances.TokenMax
                    T <: Distances.PreMetric || continue
                end
                modified_dist = mod(dist)
                @test get_max_dist(modified_dist) == d
                @test update_max_dist(modified_dist, 2*d) == mod(T(2*d))
            end
        end

    end


    @testset "`findnearest_partial` and `findall_partial` correctness with DamerauLevenshtein" begin
        ## Equal length cases

        # `d` replaced by `x`; 1 away
        str1 = "abcd"
        str2 = "abcx"
        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        @test d == Partial(DamerauLevenshtein(1))(str1, str2)
        @test matches == [(1, 1:4)] == [(d, inds)]

        # `cd` replaced by `xy`; 2 away
        str1 = "abcd"
        str2 = "abxy"

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        @test d == 2 # `max_dist + 1`
        @test isempty(matches)

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        @test d == Partial(DamerauLevenshtein(1))(str1, str2)
        @test matches == [(2, 1:4)] == [(d, inds)]

        ## Nonequal length cases

        # `d` replaced by `x`; 1 away
        str1 = "abcdef"
        str2 = "1234abcxef1234"

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        @test d == Partial(DamerauLevenshtein(1))(str1, str2)
        @test matches == [(1, 5:10)] == [(d, inds)]

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        @test d == Partial(DamerauLevenshtein(2))(str1, str2)
        @test matches == [(1, 5:10)] == [(d, inds)]

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        @test d == Partial(DamerauLevenshtein(3))(str1, str2)
        @test d == 1
        @test inds == 5:10
        @test matches == [(3, 4:9), (1, 5:10), (3, 6:11)]

        # `cde` replaced by `xyz`; 3 away
        str1 = "abcdef"
        str2 = "1234abxyzf1234"

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(1)))
        @test d == Partial(DamerauLevenshtein(1))(str1, str2)
        @test d == 2 # max_dist + 1
        @test inds == 1:0
        @test isempty(matches)

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(2)))
        @test d == Partial(DamerauLevenshtein(2))(str1, str2)
        @test d == 3 # max_dist + 1
        @test inds == 1:0
        @test isempty(matches)

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        @test d == Partial(DamerauLevenshtein(3))(str1, str2)
        @test matches == [(3, 5:10)] == [(d, inds)]

        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(4)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(4)))
        @test d == Partial(DamerauLevenshtein(4))(str1, str2)
        @test matches == [(3, 5:10)] == [(d, inds)]

        # In the first case, `cde` replaced by `xyz` (3 away); in the second, only `e` is replaced by `x` (one away)
        str1 = "abcdef"
        str2 = "1234abxyzf1234abcdxf123"
        for max_dist in (1, 2)
            d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(max_dist)))
            matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(max_dist)))
            @test d == Partial(DamerauLevenshtein(max_dist))(str1, str2)
            @test matches == [(1, 15:20)] == [(d, inds)]
        end
        # Now at 3, we find the other match.
        # We also match at "4abcdx": delete '4', substitute 'x' => 'e', and insert 'f' at the end
        # as well as "bcdxf1": insert 'a' at the start, substitute 'x' => 'e', and delete '1' at the end
        d, inds = findnearest_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        matches = findall_partial(str1, str2, Partial(DamerauLevenshtein(3)))
        @test d == Partial(DamerauLevenshtein(3))(str1, str2)
        @test (d, inds) == (1, 15:20)
        @test matches == [(3, 5:10), (3, 14:19), (1, 15:20), (3, 16:21)]
    end

    @testset "`findnearest_partial` and `findall_partial`: test other distances" begin
        # `d` replaced by `x`; 1 away
        str1 = "abcdef"
        str2 = "1234abcxef1234"

        # `Partial` unwrapping
        for d in (1,2,3), T in (DamerauLevenshtein, Levenshtein, Hamming), f in (findnearest_partial, findall_partial)
            results1 = f(str1, str2, T(d))
            results2 = f(str1, str2, Partial(T(d)))
            @test results1 == results2
        end

        @testset "$dist" for dist in (DamerauLevenshtein(), Levenshtein(), Hamming(), Jaro(), JaroWinkler(), RatcliffObershelp())
            d, inds = findnearest_partial(str1, str2, dist)
            @test inds == 5:10
            matches = findall_partial(str1, str2, dist; max_dist = d)
            @test matches == [(d, inds)]

            d1, inds1 = findnearest_partial(str1, str2, Partial(dist))
            @test d ≈ d1
            @test inds == inds1

            matches1 = findall_partial(str1, str2, Partial(dist); max_dist = d)
            @test matches1 == [(d1, inds1)]

            d2, inds2 = findnearest_partial(str1, str2, StringDistances.Normalized(dist, 1.0))
            if dist ∈ (Jaro(), JaroWinkler(), RatcliffObershelp())
                @test d2 ≈ d

                matches2 = findall_partial(str1, str2, StringDistances.Normalized(dist, 1.0); max_dist = d)
                @test matches2 == [(d2, inds2)]
                matches3 = findall_partial(str1, str2, StringDistances.Normalized(dist, d))
                @test matches3 == matches2
            else
                @test d2 * length(str1) ≈ d

                matches2 = findall_partial(str1, str2, StringDistances.Normalized(dist, 1.0); max_dist = d/length(str1) + eps())
                @test matches2 == [(d2, inds2)]
                matches3 = findall_partial(str1, str2, StringDistances.Normalized(dist, d/length(str1) + eps()))
                @test matches3 == matches2
            end
            @test inds == inds2
        end
    end

    @testset "Non-string tests" begin
        v = [6,4,1,3,2,6]
        @test findnearest_partial(1:3, v, DamerauLevenshtein()) == (1,3:5)

        matches = [(2, 2:4), (2, 3:5), (2, 4:6)]
        @test findnearest_partial(1:3, v, Hamming()) ∈ matches
        @test findall_partial(1:3, v, Hamming(), max_dist = 2) == matches
    end
end
