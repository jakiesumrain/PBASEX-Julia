using Random: Xoshiro, randexp

# Exact Poisson draws by successive arrivals of a unit-rate process.
function poisson_draw(rng, mean_count)
    arrivals = 0
    time = randexp(rng)
    while time <= mean_count
        arrivals += 1
        time += randexp(rng)
    end
    return arrivals
end

@testset "Raw-count spectrum and peak uncertainties" begin
    mktempdir() do dir
        rng = Xoshiro(87195)
        for geometry in (:quadrant, :half_plane)
            x = collect(0:5)
            orders = geometry === :quadrant ? [0, 2] : [0, 1, 2]
            basis = joinpath(dir, string(geometry) * ".h5")
            data = BASISData{Float64}(basis; x, k=[1.5, 4.0], l=orders, params=[1.4],
                                      geometry, make_images=false)
            coefficients = geometry === :quadrant ? [60.0 8; 75 -5] : [60.0 10 8; 75 -8 -5]
            folded = reshape(data.U * (data.S .* (data.V' * vec(coefficients))), :, length(x))
            mean_image = if geometry === :quadrant
                unfold_quadrant(folded; quadrant_count=1)
            else
                # Trim one side to exercise folded positions with only one observation.
                folded[:, [6:-1:2; 1:5]]
            end
            options = (; basis, geometry, center=(6, 6), alpha=0.5, reg=0.002,
                       quadrant_filter=[1, 1, 1, 0], make_images=false)
            window = (1.25, 6.25)
            # Independent trapezoids at the window endpoints and interior samples.
            function peak_components(fit)
                profiles = (data.frk' * fit.coefficients) .* x # 2 * alpha = 1
                left = (profiles[2, :] + profiles[3, :]) / 2   # E = 1.25
                right = (profiles[4, :] + profiles[5, :]) / 2  # E = 6.25
                return (2-1.25)/2 * (left + profiles[3, :]) +
                       (4.5-2)/2 * (profiles[3, :] + profiles[4, :]) +
                       (6.25-4.5)/2 * (profiles[4, :] + right)
            end
            target_fit = pbasex(mean_image; options...)
            components = peak_components(target_fit)
            target = [target_fit.IE[2:end]; components[1]; components[2:end] ./ components[1];
                      vec(target_fit.betas[2:end, :])]
            replicates = 1000
            estimates = zeros(length(target), replicates)
            errors = similar(estimates)
            for trial in 1:replicates
                counts = poisson_draw.(Ref(rng), mean_image)
                analysis = pbasex_counts(counts; options..., peak_windows=[window])
                peak = only(analysis.peaks)
                estimates[:, trial] = [analysis.result.IE[2:end]; peak.yield; peak.betas;
                                       vec(analysis.result.betas[2:end, :])]
                errors[:, trial] = [analysis.spectrum_se[2:end]; peak.yield_se; peak.beta_se;
                                    vec(analysis.beta_se[2:end, :])]
                if trial == 1
                    ordinary = pbasex(Float64.(counts); options...)
                    @test ordinary.IE ≈ analysis.result.IE
                    expected = peak_components(ordinary)
                    @test peak.yield ≈ expected[1]
                    @test peak.betas ≈ expected[2:end] ./ expected[1]
                    @test iszero(analysis.spectrum_se[1])
                    without_windows = pbasex_counts(counts; options...)
                    @test isequal(without_windows.beta_se, analysis.beta_se)
                    @test size(analysis.beta_se) == size(ordinary.betas)
                    @test all(isnan, analysis.beta_se[1, :])
                end
            end
            average = sum(estimates; dims=2) ./ replicates
            empirical_sd = sqrt.(vec(sum(abs2, estimates .- average; dims=2)) ./ (replicates-1))
            rms_se = sqrt.(vec(sum(abs2, errors; dims=2)) ./ replicates)
            ratio = rms_se ./ empirical_sd
            coverage = vec(sum(abs.(estimates .- target) .<= errors; dims=2)) ./ replicates
            @info "Poisson validation" geometry replicates se_to_sd=extrema(ratio) coverage=extrema(coverage)
            @test all(r -> 0.90 < r < 1.10, ratio)
            @test all(p -> 0.62 < p < 0.74, coverage)

            empty_image = zeros(Int, size(mean_image))
            dense_windows = [(i/2, (i+1)/2) for i in 0:10]
            hist = round.(Int, mean_image)
            dense = pbasex_counts(hist; options..., windows=dense_windows, native_errors=false)
            @test isnothing(dense.spectrum_se)
            @test isnothing(dense.beta_se)
            for (limits, batched) in zip(dense_windows, dense.windows)
                single = only(pbasex_counts(hist; options..., windows=[limits]).windows)
                @test batched.window == limits
                @test batched.yield ≈ single.yield
                @test batched.yield_se ≈ single.yield_se
                @test all(isapprox.(batched.betas, single.betas; nans=true))
                @test all(isapprox.(batched.beta_se, single.beta_se; nans=true))
            end
            zero_fit = pbasex_counts(empty_image; options..., peak_windows=[window, (0, 2)])
            @test all(iszero, zero_fit.spectrum_se)
            @test all(isnan, zero_fit.beta_se)
            @test length(zero_fit.peaks) == 2
            @test all(isnan, zero_fit.peaks[1].beta_se)
            @test isempty(pbasex_counts(empty_image; options...).peaks)
            @test_throws ArgumentError pbasex_counts(mean_image; options...)
            @test_throws ArgumentError pbasex_counts(empty_image; options..., peak_windows=[(-1, 2)])
            # An image with one occupied pixel has a single noise source.
            # Its propagated spectrum SE equals |IE| / sqrt(n), even when
            # the quadrant fold reads that same pixel more than once.
            for location in ((6, 6), (4, 6), (6, 4), (4, 5), (4, 1))
                isolated = copy(empty_image)
                isolated[location...] = 25
                single = pbasex_counts(isolated; options..., windows=dense_windows)
                @test single.spectrum_se ≈ abs.(single.result.IE) ./ 5
                # With one source of noise its scale cancels from every ratio.
                reliable = single.result.IE .> single.spectrum_se
                @test all(abs.(single.beta_se[reliable, :]) .< 1e-10)
                @test all(isnan, single.beta_se[.!reliable, :])
                for peak in single.windows
                    @test peak.yield_se ≈ abs(peak.yield) / 5 atol=1e-12
                    if peak.yield > peak.yield_se
                        @test all(abs.(peak.beta_se) .< 1e-10)
                    end
                end
            end
        end
    end
end
