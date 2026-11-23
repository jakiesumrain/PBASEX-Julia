# Unit tests for PBASEX-Julia.
#
# Run with:
#     julia --project=. -e 'using Pkg; Pkg.test()'
#
# Instantiate the project first:
#     julia --project=. -e 'using Pkg; Pkg.instantiate()'
#
using pBASEX
using pBASEX.MathUtils: meshgrid
using Test

@testset "PBASEX.jl" begin

    @testset "Legendre polynomials" begin
        # Orders 0-4 use closed forms; higher orders use Bonnet's recurrence.
        @test leg(0, 0.5) == 1.0
        @test leg(1, 0.5) == 0.5
        @test leg(2, 0.5) == -0.125          # (3x^2 - 1)/2
        @test leg(3, 0.5) == -0.4375         # (5x^3 - 3x)/2
        @test leg(4, 0.5) == -0.2890625      # (35x^4 - 30x^2 + 3)/8

        # Array method agrees with the scalar method, element-wise.
        grid = [1.0, 0.5, 0.0, -0.5, -1.0]
        for l in 0:8
            @test leg(l, grid) ≈ [leg(l, v) for v in grid]
        end

        # Known values.
        @test leg(2, [1.0, 0.0, -1.0]) ≈ [1.0, -0.5, 1.0]
        @test leg(6, 0.0) ≈ -0.3125          # P_6(0) = -5/16

        # P_l(1) == 1 and P_l(-1) == (-1)^l for every order.
        for l in 0:10
            @test leg(l, 1.0) ≈ 1.0
            @test leg(l, -1.0) ≈ (iseven(l) ? 1.0 : -1.0)
        end

        # Odd orders vanish at the origin, even ones do not.
        @test leg(5, 0.0) ≈ 0.0 atol = 1e-14
        @test leg(7, 0.0) ≈ 0.0 atol = 1e-14

        @test_throws AssertionError leg(-1, 0.5)
    end

    @testset "meshgrid" begin
        X, Y = meshgrid([1, 2, 3], [4, 5])
        @test size(X) == (2, 3)
        @test size(Y) == (2, 3)
        @test X == [1 2 3; 1 2 3]
        @test Y == [4 4 4; 5 5 5]
    end

    @testset "foldQuadrant" begin
        img = ones(5, 5)

        # All four quadrants are summed onto the first one.
        folded = foldQuadrant(img, 3, 3)
        @test size(folded) == (3, 3)
        @test all(folded .== 4.0)

        # A single selected quadrant is counted once.
        single_quadrant = foldQuadrant(img, 3, 3; quadrant_filter = [1, 0, 0, 0])
        @test size(single_quadrant) == (3, 3)
        @test all(single_quadrant .== 1.0)

        # Values add up as advertised for a non-uniform image.
        # Index [1, 1] is the folding centre, so it appears in all four
        # quadrants; index [3, 3] is the far corner, so it gathers the four
        # corners of the original image.
        graded = reshape(1.0:25.0, 5, 5)
        summed = foldQuadrant(graded, 3, 3)
        @test summed[1, 1] == 4 * graded[3, 3]
        @test summed[3, 3] ==
              graded[1, 1] + graded[1, 5] + graded[5, 1] + graded[5, 5]

        @test_throws AssertionError foldQuadrant(img, 3, 3; quadrant_filter = [1, 0])
        @test_throws AssertionError foldQuadrant(img, 0, 3)
        @test_throws AssertionError foldQuadrant(img, 3, 6)
    end

    @testset "unfoldQuadrant" begin
        folded = fill(4.0, 3, 3)

        unfolded = unfoldQuadrant(folded)          # center_offset = true
        @test size(unfolded) == (5, 5)
        @test all(unfolded .≈ 1.0)                 # divided by the 4 quadrants

        # Folding then unfolding a fourfold-symmetric image is the identity.
        original = ones(5, 5)
        @test unfoldQuadrant(foldQuadrant(original, 3, 3)) ≈ original

        @test_throws AssertionError unfoldQuadrant(ones(1, 1))
    end

    @testset "resizeFolded" begin
        grid = reshape(1.0:16.0, 4, 4)

        cropped = resizeFolded(grid, 2)
        @test size(cropped) == (2, 2)
        @test cropped == grid[1:2, 1:2]            # anchored top-left

        padded = resizeFolded(ones(2, 2), 4)
        @test size(padded) == (4, 4)
        @test padded[1:2, 1:2] == ones(2, 2)
        @test all(padded[3:4, :] .== 0.0)
        @test all(padded[:, 3:4] .== 0.0)

        @test_throws AssertionError resizeFolded(ones(3, 4), 2)
        @test_throws AssertionError resizeFolded(ones(2, 2), 0)
    end

    @testset "basis cache filenames" begin
        x = collect(0:7)
        k = [1.5, 5.5]
        l = [0, 2]

        name = generate_hashed_filename(x, k, l)
        @test startswith(name, "BASISData_")
        @test endswith(name, ".h5")
        @test length(name) == length("BASISData_") + 64 + length(".h5")  # SHA-256

        # Deterministic for identical input.
        @test name == generate_hashed_filename(x, k, l)

        # Sensitive to the grid, radial centers, and angular orders.
        @test name != generate_hashed_filename(x, k, [0, 2, 4])
        @test name != generate_hashed_filename(x, [1.5, 6.5], l)
        @test name != generate_hashed_filename(collect(0:8), k, l)

        # Prefix and suffix are honoured.
        @test startswith(generate_hashed_filename(x, k, l; prefix = "b_"), "b_")
        @test endswith(generate_hashed_filename(x, k, l; suffix = ".dat"), ".dat")
    end

    @testset "argument validation" begin
        # These all fail before any basis is computed or cached.
        @test_throws AssertionError pbasex(rand(4, 3))                       # not square
        @test_throws AssertionError pbasex(rand(4, 4); reg = -1.0)
        @test_throws AssertionError pbasex(rand(4, 4); quadrant_filter = [1, 1])
        @test_throws AssertionError pbasex(rand(4, 4); quadrant_filter = [1, 0, 0, 2])
        @test_throws AssertionError pbasex(rand(4, 4); alpha = 0.0)
        @test_throws AssertionError pbasex(rand(4, 4); center = (0, 2))
        @test_throws AssertionError pbasex(rand(4, 4); weights = ones(4, 3))
    end

    @testset "package inversion and cache round trip" begin
        mktempdir() do dir
            x = collect(0:7)
            k = [1.5, 5.5]
            l = [0, 2]
            basis = joinpath(dir, "basis.h5")
            data = BASISData{Float64}(basis; x, k, l)
            coefficients = [2.0 0.5; 1.0 -0.2]
            projection = data.U * (data.S .* (data.V' * vec(coefficients)))
            image = unfoldQuadrant(reshape(projection, length(x), length(x)))

            @test isfile(basis)
            result = pbasex(image; basis, x, k, l)
            @test result isa PBASEXResult{Float64}
            @test result.coefficients ≈ coefficients
            @test result.reconstructed_image ≈ image
            @test size(result.inverted_image) == size(image)
        end
    end

end
