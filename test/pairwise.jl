using StringDistances, Unicode, Test, Random

@testset "pairwise" begin

TestStrings = ["", "abc", "bc", "kitten"]

@testset "pairwise" begin
	for DT in [Levenshtein, Jaro]
		d = DT()
		R = pairwise(d, TestStrings)

		@test R isa Matrix{Float64, 2}
		@test size(R) == (4, 4)

		# No distance on the diagonal, since comparing strings to themselves
		@test R[1, 1] == 0.0
		@test R[2, 2] == 0.0
		@test R[3, 3] == 0.0
		@test R[4, 4] == 0.0

		# First row is comparing "" to the other strings, so:
		@test R[1, 2] == evaluate(d, "", "abc")
		@test R[1, 3] == evaluate(d, "", "bc")
		@test R[1, 4] == evaluate(d, "", "kitten")

		# Second row is comparing "abc" to the other strings, so:
		@test R[2, 3] == evaluate(d, "abc", "bc")
		@test R[2, 4] == evaluate(d, "abc", "kitten")

		# Third row row is comparing "bc" to the other strings, so:
		@test R[3, 4] == evaluate(d, "bc", "kitten")

		# Matrix is symmetric
		for i in 1:4
			for j in (i+1):4
				@test R[i, j] == R[j, i]
			end
		end
	end
end

end
