using Random: Xoshiro
using Statistics: std, mean

@testset "Detector rotation in degrees" begin
    pixels = reshape(1:24, 4, 6)
    center = (5, 2)
    for angle in (-90, 0, 90.0, 180, 270, 360, 450)
        turns = PBASEX.ImageProcessing.rotation_turns(angle)
        oriented = rotl90(pixels, round(Int, angle / 90))
        grid = PBASEX.ImageProcessing.rotated_grid(size(pixels), center, turns)
        @test size(oriented) == grid.shape
        @test oriented[reverse(grid.center)...] == pixels[reverse(center)...]
    end
end

@testset "Whole-shot bootstrap sampling uncertainty" begin
    rng = Xoshiro(915)
    shots = 120
    draws = 300
    # Two independent detector regions with equal Poisson rates. Their
    # sum is the yield and their signed difference is an angular component.
    for correlated in (false,true)
        estimates = Float64[]
        predicted = Float64[]
        covered = Bool[]
        poisson_variances = Float64[]
        for experiment in 1:160
            left = [poisson_draw(rng,6) for _ in 1:shots]
            right = [poisson_draw(rng,6) for _ in 1:shots]
            yields = left + right
            angular = Float64.(left - right)
            if correlated
                # A shared angular offset persists over three successive shots.
                angular .+= yields .* repeat(0.25 .* randn(rng,shots÷3);inner=3)
            end
            components = permutedims(hcat(yields,angular))
            pooled = vec(sum(components;dims=2))
            samples = PBASEX.ShotUncertainties.resample_totals(components,[1:shots],
                                        correlated ? 8 : 1,draws,rand(rng,UInt64))
            summary = PBASEX.ShotUncertainties.window_summary((1,2),pooled,samples,0.9,false)
            push!(estimates,only(summary.betas))
            push!(predicted,only(summary.beta_se))
            lo,hi = only(summary.beta_ci)
            push!(covered,lo<=0<=hi)
            push!(poisson_variances,(1-only(summary.betas)^2)/pooled[1])
        end
        rms_se = sqrt(mean(abs2,predicted))
        @test 0.78 < rms_se/std(estimates) < 1.20
        @test 0.78 < mean(covered) < 0.97
        if !correlated
            @test 0.95 < rms_se/sqrt(mean(poisson_variances)) < 1.05
        end
    end
end

@testset "Streamed shot estimates match direct bootstrap refits" begin
    mktempdir() do dir
        path = joinpath(dir,"shots.h5")
        ids = collect(0:2:30)
        rng = Xoshiro(31)
        # Deliberately interleave frames and leave the last frame empty.
        frames = rand(rng,ids[1:end-1],1300)
        xs,ys = 9 .* rand(rng,1300),9 .* rand(rng,1300)
        xs[1] = 9.0
        ys[1] = 0.0
        # Coordinates and alpha in a checkpoint use physical detector pixels.
        # This compact on-disk integer would overflow if squared before conversion.
        resolution = UInt8(16)
        h5open(path,"w") do f
            attributes(f)["schema_version"] = 1
            shots = create_group(f,"shots")
            shots["id"] = ids
            events = create_group(f,"events")
            events["shot_id"],events["x_corrected"],events["y_corrected"] = frames,xs ./ resolution,ys ./ resolution
            g = create_group(f,"geometry")
            g["x_edges"] = g["y_edges"] = collect(0.:9.) ./ resolution
            for (key,value) in (("center_xy",[3.5,4.5] ./ resolution),("resolution_factor",resolution),
                                ("x_direction","right"),("y_direction","down"),("coordinate_units","detector_pixel"))
                attributes(g)[key] = value
            end
        end
        weights = range(0.8,1.2;length=length(ids))
        counts = zeros(Int,9,9)
        images = [zeros(9,9) for _ in ids]
        for (id,x,y) in zip(frames,xs,ys)
            s = findfirst(==(id),ids)
            images[s][min(floor(Int,y)+1,9),min(floor(Int,x)+1,9)] += weights[s]
            counts[min(floor(Int,y)+1,9),min(floor(Int,x)+1,9)] += 1
        end
        images = rotl90.(images)
        windows = [(1.,9.),(4.,16.),(1.,9.)]
        for geometry in (:quadrant,:half_plane)
            basis = joinpath(dir,string(geometry)*".h5")
            options = (; basis,geometry,center=(5,6),quadrant_filter=[1,1,1,0],
                       reg=0.1,make_images=false)
            orders = geometry === :quadrant ? [0,2] : [0,1,2]
            pbasex(sum(images);options...,x=0:4,k=[1.,3.],l=orders)
            data = BASISData(basis;make_images=false)
            function components(image,w)
                fit = pbasex(image;options...)
                angular = (data.frk' * fit.coefficients) .* (data.x ./ 2)
                # Endpoints coincide with native samples in this fixture.
                idx = findall(e -> w[1]<=e<=w[2],fit.E)
                return vec(sum((angular[idx[1:end-1],:] + angular[idx[2:end],:]) .*
                               (diff(fit.E[idx]) ./ 2);dims=1))
            end
            settings = (; basis,windows,rotation=90,alpha=Float64(resolution)^2,quadrant_filter=options.quadrant_filter,
                        reg=options.reg,shot_weights=weights,keep_replicates=true,replicates=12,
                        block_length=2,run_ids=repeat([1,2];inner=8))
            result = pbasex_shots(path;settings...,rng=Xoshiro(99),event_chunk=57,
                                  window_batch=2,basis_chunk_bytes=512)
            for reuse_shot_work in (false, true)
                resident = pbasex_shots(path;settings...,rng=Xoshiro(99),event_chunk=57,
                                        window_batch=2,resident_u=true,reuse_shot_work)
                for (streamed, loaded) in zip(result.windows, resident.windows)
                    @test loaded.replicates.yield ≈ streamed.replicates.yield
                    @test all(isapprox.(loaded.replicates.betas, streamed.replicates.betas; nans=true))
                    @test loaded.valid_replicates == streamed.valid_replicates
                end
            end
            @test result.shots == 16
            legacy = pbasex_shots(path; settings..., rotation=nothing, quarter_turns=1,
                                  rng=Xoshiro(99), event_chunk=57, window_batch=2, basis_chunk_bytes=512)
            @test legacy.windows[1].yield == result.windows[1].yield
            @test_throws ArgumentError pbasex_shots(path; settings..., quarter_turns=1)
            @test_throws ArgumentError pbasex_shots(path; settings..., rotation=45)
            @test_throws ArgumentError pbasex_counts(counts; options..., rotation=45)
            pooled = pbasex_counts(counts; options..., windows, center=(4,5), rotation=90)
            explicit = pbasex_counts(rotl90(counts); options..., windows)
            @test pooled.result.IE ≈ explicit.result.IE
            @test pooled.spectrum_se ≈ explicit.spectrum_se
            shots = pbasex_shots(path; settings..., geometry, center=(4,5), shot_weights=nothing, rng=Xoshiro(99))
            @test pooled.legendre_orders == shots.legendre_orders
            @test pooled.peaks === pooled.windows
            @test length(pooled.windows) == length(shots.windows) == length(windows)
            for (a,b) in zip(pooled.windows,shots.windows)
                @test a.window == b.window
                @test a.yield ≈ b.yield
                @test a.spectrum ≈ b.spectrum
                @test a.betas ≈ b.betas
                @test size(a.beta_se) == size(b.beta_se)
                for bin in (a,b)
                    @test bin.spectrum_se ≈ bin.yield_se / (bin.window[2]-bin.window[1])
                end
            end
            @test_throws ArgumentError pbasex_counts(counts; options..., windows, peak_windows=windows)
            @test_throws ArgumentError pbasex_shots(path; settings..., geometry=:mismatch)
            @test_throws ArgumentError pbasex_shots(path;settings...,block_length=8)
            @test result.legendre_orders == orders
            boot = Xoshiro(rand(Xoshiro(99),UInt64))
            for trial in 0:12
                draw = collect(1:16)
                if trial > 0
                    empty!(draw)
                    for first in (1,9), _ in 1:4
                        start = rand(boot,first:first+6)
                        append!(draw,start:start+1)
                    end
                end
                image = sum(images[draw])
                for (w,summary) in zip(windows,result.windows)
                    b = components(image,w)
                    yield = trial==0 ? summary.yield : summary.replicates.yield[trial]
                    betas = trial==0 ? summary.betas : summary.replicates.betas[:,trial]
                    @test yield ≈ b[1]
                    @test betas ≈ b[2:end] ./ b[1]
                end
            end
            @test result.windows[1].replicates.betas ≈ result.windows[3].replicates.betas
            for reuse_shot_work in (false, true)
                zero_result = pbasex_shots(path;settings...,shot_weights=zeros(16),
                                           resident_u=reuse_shot_work,reuse_shot_work)
                @test all(w -> w.valid_replicates == 0 && all(isnan, w.betas) &&
                               all(isnan, w.beta_se) && all(isnan, w.replicates.betas),
                          zero_result.windows)
            end
            streamed = []
            second = pbasex_shots(path;settings...,rng=Xoshiro(99),keep_replicates=false,
                                  replicate_sink=(i,draws)->push!(streamed,draws))
            @test all(w -> isnothing(w.replicates),second.windows)
            @test streamed[2].yield ≈ result.windows[2].replicates.yield
            if geometry === :quadrant
                matlab_basis = joinpath(dir,"matlab.h5")
                h5open(matlab_basis,"w") do f
                    for name in (:x,:k,:l,:S,:frk)
                        f[string(name)] = getproperty(data,name)
                    end
                    f["Up"] = data.U
                    f["V"] = permutedims(data.V)
                end
                matlab = pbasex_shots(path;settings...,basis=matlab_basis,rng=Xoshiro(99))
                @test matlab.windows[2].replicates.betas ≈ result.windows[2].replicates.betas
                h5open(matlab_basis,"r+") do f
                    delete_object(f,"frk")
                    f["frk"] = data.frk[:,1:1]
                end
                @test_throws ArgumentError pbasex_shots(path; settings..., basis=matlab_basis)
            end
        end
    end
end
