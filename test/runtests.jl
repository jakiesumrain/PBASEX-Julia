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
    include("counting_errors.jl")
    include("shot_uncertainties.jl")

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

    @testset "argument validation" begin
        @test_throws UndefKeywordError pbasex(ones(5, 5))
        mktempdir() do dir
            basis = joinpath(dir, "unused.h5")
            @test_throws ArgumentError pbasex(ones(5, 5); basis=dir)
            @test_throws ArgumentError BASISData(dir)
            @test_throws ArgumentError BASISData(joinpath(dir, "missing", "basis.h5"))
            @test_throws AssertionError pbasex(rand(4, 3); basis)
            @test_throws AssertionError pbasex(rand(4, 4); basis, reg=-1.0)
            @test_throws AssertionError pbasex(rand(4, 4); basis, quadrant_filter=[1, 1])
            @test_throws AssertionError pbasex(rand(4, 4); basis, quadrant_filter=[1, 0, 0, 2])
            @test_throws AssertionError pbasex(rand(4, 4); basis, alpha=0.0)
            @test_throws AssertionError pbasex(rand(4, 4); basis, center=(0, 2))
            @test_throws AssertionError pbasex(rand(4, 4); basis, weights=ones(4, 3))
            @test_throws ArgumentError pbasex(ones(5, 5); basis, l=[2, 4])
            @test_throws ArgumentError pbasex(ones(5, 5); basis, l=[0, 1])
            @test_throws ArgumentError pbasex(ones(5, 5); basis, l=[0, 0])
            @test isempty(readdir(dir))
        end
    end

    @testset "package inversion and basis file round trip" begin
        mktempdir() do dir
            x = collect(0:7)
            k = [1.5, 5.5]
            l = [0, 2]
            basis = joinpath(dir, "basis.h5")
            data = BASISData{Float64}(basis; x, k, l)
            coefficients = [2.0 0.5; 1.0 -0.2]
            projection = data.U * (data.S .* (data.V' * vec(coefficients)))
            image = unfold_quadrant(reshape(projection, length(x), length(x)); quadrant_count=1)

            @test isfile(basis)
            result = pbasex(image; basis)
            @test result isa PBASEXResult{Float64}
            @test result.coefficients ≈ coefficients
            half = pbasex(image; basis=joinpath(dir, "half.h5"), x, k, l, geometry=:half_plane)
            @test half.coefficients ≈ result.coefficients
            @test half.IE ≈ result.IE
            @test half.betas[2:end, :] ≈ result.betas[2:end, :]
            @test half.inverted_image ≈ result.inverted_image
            @test half.reconstructed_image ≈ result.reconstructed_image
            @test maximum(abs, result.residuals) < 1e-12
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

    @testset "explicit basis files and lazy slice" begin
        mktempdir() do dir
            x, l = 0:7, [0, 2]
            basis = joinpath(dir, "gaussian-width-1p4.h5")
            in_memory = BASISData(; x, l, params=[1.4], make_images=false)
            @test isempty(readdir(dir))
            first_basis = BASISData(basis; x, l, params=[1.4], make_images=false)
            @test first_basis.Ginv === nothing
            @test first_basis.S == in_memory.S
            @test readdir(dir) == [basename(basis)]
            @test !h5open(f -> haskey(f, "Ginv"), basis, "r")

            @test_throws ArgumentError BASISData(basis; params=[3.0], make_images=false)
            other_width = BASISData(joinpath(dir, "gaussian-width-3.h5");
                                    x, l, params=[3.0], make_images=false)
            @test other_width.frk != first_basis.frk
            weighted = BASISData(joinpath(dir, "gaussian-weighted.h5");
                                 x, l, params=[1.4], weights=fill(2.0, 8, 8), make_images=false)
            @test weighted.S ≈ 2 .* first_basis.S

            complete = BASISData(basis)
            @test complete.S == first_basis.S
            @test size(complete.Ginv) == (64, 4)
            @test !h5open(f -> haskey(f, "Ginv"), basis, "r")

            radial = (r, k, p) -> exp(-abs(r - k) / p[1])
            custom_path = joinpath(dir, "custom.h5")
            custom = BASISData(custom_path; x, l, rBF=radial, make_images=false)
            @test BASISData(custom_path; make_images=false).S == custom.S
            @test_throws ArgumentError BASISData(custom_path; rBF=radial)

            tagged_path = joinpath(dir, "custom-tagged.h5")
            tagged = BASISData(tagged_path; x, l, rBF=radial, cache_key="exponential-v1", make_images=false)
            reopened = BASISData(tagged_path; rBF=radial, cache_key="exponential-v1")
            @test reopened.S == tagged.S
            @test size(reopened.Ginv) == (64, 4)
            @test_throws ArgumentError BASISData(tagged_path; rBF=radial, cache_key="exponential-v2")

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

    @testset "Axial symmetry with odd Legendre orders" begin
        geometry = :half_plane
        x, k, l = 0:5, [0.0, 1.0, 2.0], [0, 1, 2]
        side = length(x)
        center = (side, side)
        width, asymmetry, amplitude = 2.0, 0.35, 1.7
        params = [width]
        # For this analytic example, the radial basis indices encode powers of r/s.
        radial(r, power, parameters) = (r / parameters[1])^Int(power) * exp(-r^2 / (2 * parameters[1]^2))
        limit(R, power, parameters) = sqrt(max(0.0, (14 * parameters[1])^2 - R^2))
        slice = [amplitude * exp(-(y^2 + z^2) / (2 * width^2)) * (1 + asymmetry * z / width)^2
                 for z in last(x):-1:-last(x), y in -last(x):last(x)]
        image = (sqrt(2π) * width) .* slice
        expected = zeros(3, 3)
        expected[1, 1] = 4π * amplitude
        expected[3, 1] = 4π * amplitude * asymmetry^2 / 3
        expected[2, 2] = 8π * amplitude * asymmetry
        expected[3, 3] = 8π * amplitude * asymmetry^2 / 3

        mktempdir() do dir
            basis = joinpath(dir, "signed-axis.h5")
            result = pbasex(image; basis, geometry, x, k, l, params, rBF=radial, zIP=limit)
            @test result.legendre_orders == l
            @test result.coefficients ≈ expected
            @test result.reconstructed_image ≈ image
            @test result.inverted_image ≈ slice
            @test maximum(abs, result.residuals) < 1e-10
            @test h5read(basis, "geometry") == "half_plane"
            @test BASISData(basis).geometry === geometry
            @test_throws ArgumentError pbasex(image; basis)
            @test_throws ArgumentError pbasex(image; basis, geometry, quadrant_filter=[1, 1, 0, 0])
            @test_throws ArgumentError BASISData(basis; geometry=:quadrant)

            q = asymmetry .* collect(x[2:end]) ./ width
            scale = 1 .+ q.^2 ./ 3
            @test result.betas[2:end, 1] ≈ 2 .* q ./ scale
            @test result.betas[2:end, 2] ≈ 2 .* q.^2 ./ (3 .* scale)
            @test all(isnan, result.betas[1, :])

            flipped = pbasex(reverse(image; dims=1); basis, geometry, make_images=false)
            @test flipped.IE ≈ result.IE
            @test flipped.betas[2:end, 1] ≈ -result.betas[2:end, 1]
            @test flipped.betas[2:end, 2] ≈ result.betas[2:end, 2]
            symmetric = pbasex((image + reverse(image; dims=1)) ./ 2; basis, geometry, make_images=false)
            @test maximum(abs, symmetric.betas[2:end, 1]) < 1e-10

            # Retain the original pixel coordinates after cropping and ignore masked data.
            cropped = image[2:end, 1:end-2]
            weights = ones(size(cropped))
            weights[1, 5] = 0
            cropped[1, 5] = NaN
            weighted_basis = joinpath(dir, "signed-weighted.h5")
            weighted = pbasex(cropped; basis=weighted_basis, geometry, x, k, l, params,
                              rBF=radial, zIP=limit, weights, center=(center[1], center[2]-1),
                              make_images=false)
            @test weighted.coefficients ≈ expected
            @test any(iszero, BASISData(weighted_basis; make_images=false).weights)

            legacy = joinpath(dir, "quadrant.h5")
            BASISData(legacy; x, make_images=false)
            @test_throws ArgumentError pbasex(image; basis=legacy, geometry)
        end

        for radial_basis in (:gaussian, :cubic_bspline)
            data = BASISData(; geometry, x, k=[1.5, 3.5], l, radial_basis, params=[1.0])
            G = data.U * (data.S .* data.V')
            for basis_matrix in (G, data.Ginv)
                even = reshape(basis_matrix[:, 1], 2 * side - 1, side)
                odd = reshape(basis_matrix[:, 3], 2 * side - 1, side)
                @test reverse(even; dims=1) ≈ even
                @test reverse(odd; dims=1) ≈ -odd
                @test maximum(abs, odd[side, :]) < 1e-13
            end
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
            BASISData(joinpath(dir, "gaussian.h5"); x, k, l, params, make_images=false)
            spline_path = joinpath(dir, "spline.h5")
            cached = BASISData(spline_path; x, k, l, params, radial_basis=family, make_images=false)
            @test sort(readdir(dir)) == ["gaussian.h5", "spline.h5"]
            @test cached.Ginv === nothing
            @test cached.frk == spline.frk
            restored = BASISData(spline_path; x, k, l, params, radial_basis=family)
            @test restored.S == cached.S
            @test restored.Ginv == spline.Ginv

            c = [2.0 0.5; 1.0 -0.2]
            projection = spline.U * (spline.S .* (spline.V' * vec(c)))
            image = unfold_quadrant(reshape(projection, length(x), length(x)); quadrant_count=1)
            result = pbasex(image; basis=spline_path, x, k, l, params, radial_basis=family,
                            make_images=false)
            @test result.coefficients ≈ c
            @test pbasex(image; basis=spline_path).inverted_image ≈
                  unfold_quadrant(reshape(spline.Ginv * vec(c), length(x), length(x)); quadrant_count=1) ./ (4π)
            @test_throws ArgumentError pbasex(image; basis=spline_path, radial_basis=:gaussian)
        end
    end

end
