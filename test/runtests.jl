# Unit tests for PBASEX-Julia.
#
# Run with:
#     julia --project=. -e 'using Pkg; Pkg.test()'
#
# Instantiate the project first:
#     julia --project=. -e 'using Pkg; Pkg.instantiate()'
#
using PBASEX
using PBASEX.MathUtils: meshgrid
using HDF5
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

    @testset "fold_quadrants" begin
        img = ones(5, 5)

        # All four quadrants are summed onto the first one.
        folded = fold_quadrants(img, 3, 3)
        @test size(folded) == (3, 3)
        @test all(folded .== 4.0)

        # A single selected quadrant is counted once.
        single_quadrant = fold_quadrants(img, 3, 3; quadrant_filter = [1, 0, 0, 0])
        @test size(single_quadrant) == (3, 3)
        @test all(single_quadrant .== 1.0)
        labeled = [11 12 13; 21 22 23; 31 32 33]
        for (i, corner) in enumerate((13, 11, 31, 33))
            filter = [Int(i == j) for j in 1:4]
            @test fold_quadrants(labeled, 2, 2; quadrant_filter=filter)[2, 2] == corner
        end
        @test_throws AssertionError fold_quadrants(img, 3, 3; quadrant_filter=zeros(Int, 4))

        # Values add up as advertised for a non-uniform image.
        # Index [1, 1] is the folding centre, so it appears in all four
        # quadrants; index [3, 3] is the far corner, so it gathers the four
        # corners of the original image.
        graded = reshape(1.0:25.0, 5, 5)
        summed = fold_quadrants(graded, 3, 3)
        @test summed[1, 1] == 4 * graded[3, 3]
        @test summed[3, 3] ==
              graded[1, 1] + graded[1, 5] + graded[5, 1] + graded[5, 5]

        @test_throws AssertionError fold_quadrants(img, 3, 3; quadrant_filter = [1, 0])
        @test_throws AssertionError fold_quadrants(img, 0, 3)
        @test_throws AssertionError fold_quadrants(img, 3, 6)
    end

    @testset "unfold_quadrant" begin
        folded = fill(4.0, 3, 3)

        unfolded = unfold_quadrant(folded)          # center_offset = true
        @test size(unfolded) == (5, 5)
        @test all(unfolded .≈ 1.0)                 # divided by the 4 quadrants

        # Folding then unfolding a fourfold-symmetric image is the identity.
        original = ones(5, 5)
        @test unfold_quadrant(fold_quadrants(original, 3, 3)) ≈ original
        @test unfold_quadrant(ones(3, 3); quadrant_count=1) == original
        @test unfold_quadrant(ones(3, 3, 2); center_offset=false, quadrant_count=1) == ones(6, 6, 2)
        @test_throws AssertionError unfold_quadrant(folded; quadrant_count=0)

        @test_throws AssertionError unfold_quadrant(ones(1, 1))
    end

    @testset "resize_folded_quadrant" begin
        grid = reshape(1.0:16.0, 4, 4)

        cropped = resize_folded_quadrant(grid, 2)
        @test size(cropped) == (2, 2)
        @test cropped == grid[1:2, 1:2]            # anchored top-left

        padded = resize_folded_quadrant(ones(2, 2), 4)
        @test size(padded) == (4, 4)
        @test padded[1:2, 1:2] == ones(2, 2)
        @test all(padded[3:4, :] .== 0.0)
        @test all(padded[:, 3:4] .== 0.0)

        @test resize_folded_quadrant(ones(2, 4), 3) == [1 1 1; 1 1 1; 0 0 0]
        @test resize_folded_quadrant(ones(2, 4, 2), 3)[:, :, 2] == [1 1 1; 1 1 1; 0 0 0]
        @test_throws AssertionError resize_folded_quadrant(ones(2, 2), 0)
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
        @test_throws ArgumentError pbasex(ones(5, 5); l=[2, 4])
        @test_throws ArgumentError pbasex(ones(5, 5); l=[0, 1])
        @test_throws ArgumentError pbasex(ones(5, 5); l=[0, 0])
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
            image = unfold_quadrant(reshape(projection, length(x), length(x)))

            @test isfile(basis)
            result = pbasex(image; basis)
            @test result isa PBASEXResult{Float64}
            @test result.coefficients ≈ coefficients
            @test result.reconstructed_image ≈ image
            @test size(result.inverted_image) == size(image)
            @test pbasex(image; basis, quadrant_filter=[1, 0, 0, 0]).IE ≈ result.IE
            @test pbasex(image; basis, quadrant_filter=[1, 0, 0, 0]).reconstructed_image ≈ image
            shifted = @test_logs (:warn, r"Residuals require") pbasex(image; basis, center=(7, 8))
            @test shifted.residuals === nothing
            @test_throws ArgumentError pbasex(image; basis, x=1:8)
            @test_throws ArgumentError pbasex(image; basis, params=[1.0])
            @test_throws ArgumentError pbasex(image; basis, weights=ones(size(image)))

            matlab_basis = joinpath(dir, "matlab.h5")
            h5open(matlab_basis, "w") do f
                write(f, "x", reshape(Float64.(x), :, 1))
                write(f, "k", reshape(k, 1, :))
                write(f, "l", reshape(Float64.(l), 1, :))
                write(f, "Up", data.U)
                write(f, "S", reshape(data.S, 1, :))
                write(f, "V", permutedims(data.V))
                write(f, "frk", data.frk)
                write(f, "Ginv", permutedims(data.Ginv))
            end
            reference = pbasex(image; basis=matlab_basis)
            @test reference.coefficients ≈ result.coefficients
            @test reference.inverted_image ≈ result.inverted_image

            legacy_basis = joinpath(dir, "legacy.h5")
            h5open(legacy_basis, "w") do f
                for name in (:x, :k, :l, :U, :S, :V, :frk, :Ginv)
                    write(f, String(name), getproperty(data, name))
                end
            end
            legacy = @test_logs (:warn, r"Legacy Julia basis") pbasex(image; basis=legacy_basis)
            @test legacy.IE ≈ result.IE
        end
    end

    @testset "cache identity and lazy slice" begin
        mktempdir() do dir
            x, l = 0:7, [0, 2]
            first_basis = BASISData(dir; x, l, params=[1.4], make_images=false)
            @test first_basis.Ginv === nothing
            @test !h5open(f -> haskey(f, "Ginv"), only(readdir(dir; join=true)), "r")
            other_width = BASISData(dir; x, l, params=[3.0], make_images=false)
            @test other_width.frk != first_basis.frk
            weighted = BASISData(dir; x, l, params=[1.4], weights=fill(2.0, 8, 8), make_images=false)
            @test weighted.S ≈ 2 .* first_basis.S
            @test length(readdir(dir)) == 3
            complete = BASISData(dir; x, l, params=[1.4])
            @test complete.S == first_basis.S
            @test size(complete.Ginv) == (64, 4)

            radial = (r, k, p) -> exp(-abs(r - k) / p[1])
            BASISData(dir; x, l, rBF=radial, make_images=false)
            @test length(readdir(dir)) == 3
            a = BASISData(dir; x, l, rBF=radial, cache_key="exponential-v1", make_images=false)
            b = BASISData(dir; x, l, rBF=radial, cache_key="exponential-v1", make_images=false)
            @test a.S == b.S
            @test length(readdir(dir)) == 4

            image = ones(15, 15)
            path = joinpath(dir, "weighted.h5")
            fitting = pbasex(image; basis=path, x, l, weights=image, make_images=false)
            reloaded = pbasex(image; basis=path, weights=image, make_images=false)
            @test fitting.coefficients ≈ reloaded.coefficients
            @test_throws ArgumentError pbasex(image; basis=path, make_images=false)
        end
    end

    @testset "Float32 and rank-deficient inversions" begin
        mktempdir() do dir
            image = ones(Float32, 15, 15)
            result = pbasex(image; basis=joinpath(dir, "single.h5"), x=0:7, l=[0, 2])
            @test result isa PBASEXResult{Float32}
            @test all(isfinite, result.coefficients)
            @test all(isfinite, result.inverted_image)

            zero_radial = (r, k, p) -> zero(r)
            result = pbasex(image; basis=joinpath(dir, "zero.h5"), x=0:7, l=[0, 2],
                            rBF=zero_radial, make_images=false)
            @test all(iszero, result.coefficients)
            @test all(isnan, result.betas)
        end
    end

    @testset "compact cubic B-spline basis" begin
        family = :cubic_bspline
        # With R=0 and l=0, an entire cardinal spline has area h along the ray.
        # This shell is narrow and far beyond the measured detector radius.
        h = 0.25
        distant = BASISData{Float64}(; x=0:3, k=[509.5], l=[0], params=[h],
                                     radial_basis=family, make_images=false)
        forward = distant.U * (distant.S .* distant.V')
        @test forward[1, 1] ≈ h / (2π)
        @test all(iszero, distant.frk)

        x, k, l, params = 0:7, [1.5, 5.5], [0, 2], [1.0]
        spline = BASISData{Float64}(; x, k, l, params, radial_basis=family)
        @test spline.rBF(1.5, 1.5, params) ≈ 2/3
        @test spline.rBF(2.5, 1.5, params) ≈ 1/6
        @test spline.rBF(3.5, 1.5, params) == 0
        @test spline.frk[1, 5] == 0
        @test spline.Ginv[5, 1] == 0
        @test spline.frk[1, 1] ≈ 1/48
        centered = BASISData{Float32}(; x=0:3, k=[0], l=[0], params=[1], radial_basis=family)
        @test eltype(centered.Ginv) === Float32
        @test (centered.U * (centered.S .* centered.V'))[1, 1] ≈ 1f0 / (4f0 * Float32(π))
        @test_throws ArgumentError BASISData{Float64}(; x, radial_basis=:unknown)
        @test_throws ArgumentError BASISData{Float64}(; x, radial_basis=family, params=[0])

        mktempdir() do dir
            BASISData(dir; x, k, l, params, make_images=false)
            cached = BASISData(dir; x, k, l, params, radial_basis=family, make_images=false)
            @test length(readdir(dir)) == 2
            @test cached.Ginv === nothing
            @test cached.frk == spline.frk
            restored = BASISData(dir; x, k, l, params, radial_basis=family)
            @test restored.S == cached.S
            @test restored.Ginv == spline.Ginv

            c = [2.0 0.5; 1.0 -0.2]
            projection = spline.U * (spline.S .* (spline.V' * vec(c)))
            image = unfold_quadrant(reshape(projection, length(x), length(x)))
            path = joinpath(dir, "spline.h5")
            result = pbasex(image; basis=path, x, k, l, params, radial_basis=family,
                            make_images=false)
            @test result.coefficients ≈ c
            @test pbasex(image; basis=path).inverted_image ≈
                  unfold_quadrant(reshape(spline.Ginv * vec(c), length(x), length(x)))
            @test_throws ArgumentError pbasex(image; basis=path, radial_basis=:gaussian)
        end
    end

end
