using StringDistances, Unicode, Test, Random
@testset "pairwise" begin

TestStrings1 = ["", "abc", "bc", "kitten"]
TestStrings2 = ["mew", "ab"]

TestStrings1missing = ["", "abc", "bc", missing]
TestStrings2missing = ["mew", missing]

@testset "pairwise" begin
	for DT in [Jaro, Levenshtein, DamerauLevenshtein, RatcliffObershelp,
				QGram, Cosine, Jaccard, SorensenDice, Overlap]

		d = (DT <: QGramDistance) ? DT(2) : DT()
		R = pairwise(d, TestStrings1)

		@test size(R) == (4, 4)

		# No distance on the diagonal, since comparing strings to themselves
		@test R[1, 1] == 0.0
		@test R[2, 2] == 0.0
		@test R[3, 3] == 0.0
		@test R[4, 4] == 0.0

		# Since the distance might be NaN:
		equalorNaN(x, y) = (x == y) || (isnan(x) && isnan(y))

		# First row is comparing "" to the other strings, so:
		@test equalorNaN(R[1, 2], evaluate(d, "", "abc"))
		@test equalorNaN(R[1, 3], evaluate(d, "", "bc"))
		@test equalorNaN(R[1, 4], evaluate(d, "", "kitten"))

		# Second row is comparing "abc" to the other strings, so:
		@test equalorNaN(R[2, 3], evaluate(d, "abc", "bc"))
		@test equalorNaN(R[2, 4], evaluate(d, "abc", "kitten"))

		# Third row row is comparing "bc" to the other strings, so:
		@test equalorNaN(R[3, 4], evaluate(d, "bc", "kitten"))

		# Matrix is symmetric
		for i in 1:4
			for j in (i+1):4
				@test equalorNaN(R[i, j], R[j, i])
			end
		end

		# Test also the assymetric version
		R2 = pairwise(d, TestStrings1, TestStrings2)
		@test size(R2) == (4, 2)

		@test equalorNaN(R2[1, 1], evaluate(d, "", "mew"))
		@test equalorNaN(R2[1, 2], evaluate(d, "", "ab"))

		@test equalorNaN(R2[2, 1], evaluate(d, "abc", "mew"))
		@test equalorNaN(R2[2, 2], evaluate(d, "abc", "ab"))

		@test equalorNaN(R2[3, 1], evaluate(d, "bc", "mew"))
		@test equalorNaN(R2[3, 2], evaluate(d, "bc", "ab"))

		@test equalorNaN(R2[4, 1], evaluate(d, "kitten", "mew"))
		@test equalorNaN(R2[4, 2], evaluate(d, "kitten", "ab"))

		R3 = pairwise(d, TestStrings2, TestStrings1)
		@test size(R3) == (2, 4)

		for i in 1:length(TestStrings1)
			for j in 1:length(TestStrings2)
				@test equalorNaN(R2[i, j], R3[j, i])
			end
		end

		# Ensure same result if preprocessing for QGramDistances
		if DT <: QGramDistance
			R4 = pairwise(d, TestStrings1; preprocess = true)
			@test typeof(R4) == typeof(R)
			@test size(R4) == size(R)
			for i in 1:size(R4, 1)
				for j in 1:size(R4, 2)
					@test equalorNaN(R4[i, j], R[i, j])
				end
			end
		end
		# ensures missing
		R5 = pairwise(d, TestStrings1missing; preprocess = true)
		@test eltype(R5) == Union{result_type(d, String, String), Missing}
	end
end

end
