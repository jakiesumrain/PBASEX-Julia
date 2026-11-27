module ShotUncertainties

using HDF5
using LinearAlgebra: mul!
using Random: AbstractRNG, default_rng, Xoshiro
using Statistics: std, quantile
using ..Inversion: inverse_filter, validate_angular_orders
using ..CountingErrors: window_integrals
using ..MathUtils: detector_dimensions, axial_coordinates
using ..ImageProcessing: fold_quadrants, rotation_turns, rotated_grid

# For each physical pixel, store its unique folded row and its contribution
# to that row. Multiple occurrences at quadrant axes add before averaging.
function fold_lookup(shape, center, side, geometry, filter)
    grid = detector_dimensions(side, geometry)
    index = zeros(Int, shape)
    fraction = zeros(Float64, shape)
    covered = falses(grid)
    cx, cy = center
    top_right, top_left, bottom_left, bottom_right = isone.(filter)
    directions = ((1,-1), (-1,-1), (-1,1), (1,1))
    selected = findall(isone, filter)
    if geometry === :quadrant
        shape[1] == shape[2] || throw(ArgumentError("Quadrant geometry requires a square image"))
        extent = size(fold_quadrants(zeros(UInt8, shape), cx, cy; quadrant_filter=filter))
        for q in selected, col in 1:min(side,extent[2]), row in 1:min(side,extent[1])
            pixel = CartesianIndex(cy + directions[q][2]*(row-1), cx + directions[q][1]*(col-1))
            index[pixel] = row + (col-1)*side
            fraction[pixel] += 1/length(selected)
        end
        return index, fraction, nothing
    end
    (top_left || top_right) && (bottom_left || bottom_right) ||
        throw(ArgumentError("Select at least one quadrant in each hemisphere"))
    for row in axes(covered,1), col in axes(covered,2)
        y = cy + row - side
        left = row == side ? top_left || bottom_left : (row < side ? top_left : bottom_left)
        right = row == side ? top_right || bottom_right : (row < side ? top_right : bottom_right)
        xl, xr = cx-col+1, cx+col-1
        left &= checkbounds(Bool, index, y, xl)
        right &= checkbounds(Bool, index, y, xr)
        left &= !(col == 1 && right)
        count = left + right
        count == 0 && continue
        covered[row,col] = true
        for (use, x) in ((left,xl),(right,xr))
            use || continue
            index[y,x] = row + (col-1)*grid[1]
            fraction[y,x] = 1/count
        end
    end
    return index, fraction, all(covered) ? nothing : covered
end

function pixel_bin(edges, x)
    first(edges) <= x <= last(edges) || throw(ArgumentError("Event or center outside histogram edges"))
    return min(searchsortedlast(edges, x), length(edges)-1)
end

function checkpoint_layout(file, turns, center)
    read_attribute(file, "schema_version") == 1 || throw(ArgumentError("Expected checkpoint schema 1"))
    geom = file["geometry"]
    read_attribute(geom,"x_direction") == "right" && read_attribute(geom,"y_direction") == "down" &&
        read_attribute(geom,"coordinate_units") == "detector_pixel" ||
        throw(ArgumentError("Expected detector pixels with x right and y down"))
    xe, ye = read(geom,"x_edges"), read(geom,"y_edges")
    resolution = Float64(read_attribute(geom,"resolution_factor"))
    isfinite(resolution) && resolution > 0 || throw(ArgumentError("Invalid histogram resolution"))
    for edges in (xe, ye)
        length(edges) >= 2 && all(d -> d ≈ inv(resolution), diff(edges)) ||
            throw(ArgumentError("Expected uniform edges in detector pixels"))
    end
    nx, ny = length(xe)-1, length(ye)-1
    if isnothing(center)
        x,y = read_attribute(geom,"center_xy")
        cx,cy = pixel_bin(xe,x),pixel_bin(ye,y)
        center = (cx,cy)
    end
    grid = rotated_grid((ny,nx), center, turns)
    return (; xe,ye,resolution,grid...)
end

function shot_totals(file, layout, ids, shot_index, weights, lookup, fraction, responses, chunk)
    events = file["events"]
    shots,x,y = (events[key] for key in ("shot_id","x_corrected","y_corrected"))
    size(shots) == size(x) == size(y) && ndims(x) == 1 && eltype(shots)<:Integer ||
        throw(ArgumentError("Expected aligned event columns and integer shot IDs"))
    totals = zeros(size(responses,2), length(ids))
    for start in 1:chunk:length(x)
        rows = start:min(length(x),start+chunk-1)
        frames,xs,ys = shots[rows],x[rows],y[rows]
        for i in eachindex(frames)
            s = get(shot_index, frames[i], 0)
            s > 0 || throw(ArgumentError("Event references unknown shot"))
            p = CartesianIndex(pixel_bin(layout.ye,ys[i]),pixel_bin(layout.xe,xs[i]))
            row = lookup[p]
            row == 0 && continue
            scale = weights[s]*fraction[p]
            @inbounds for j in axes(totals,1)
                totals[j,s] += scale*responses[row,j]
            end
        end
    end
    return totals
end

function resample_totals(totals, runs, block, draws, seed)
    rng = Xoshiro(seed)
    samples = zeros(size(totals,1),draws)
    multiplicity = zeros(size(totals,2))
    for draw in 1:draws
        fill!(multiplicity,0)
        for run in runs
            remaining = length(run)
            while remaining > 0
                start = rand(rng,first(run):last(run)-block+1)
                take = min(block,remaining)
                @views multiplicity[start:start+take-1] .+= 1
                remaining -= take
            end
        end
        mul!(view(samples,:,draw),totals,multiplicity)
    end
    return samples
end

function window_summary(window, pooled, samples, confidence, keep)
    yield_samples = samples[1,:]
    valid = yield_samples .> 0
    tail = (1-confidence)/2
    interval(v) = Tuple(quantile(v,[tail,1-tail]))
    yield_ci = interval(yield_samples)
    betas = pooled[2:end] ./ (pooled[1] > 0 ? pooled[1] : NaN)
    beta_draws = samples[2:end,:] ./ transpose(yield_samples)
    beta_draws[:,.!valid] .= NaN
    beta_se = fill(NaN,length(betas))
    beta_ci = fill((NaN,NaN),length(betas))
    if pooled[1] > 0 && yield_ci[1] > 0 && count(valid) >= 2
        for j in eachindex(betas)
            ratios = beta_draws[j,valid]
            beta_se[j] = std(ratios)
            beta_ci[j] = interval(ratios)
        end
    end
    width = window[2]-window[1]
    yield_se = std(yield_samples)
    replicates = keep ? (; yield=yield_samples,betas=beta_draws) : nothing
    return (; window,yield=pooled[1],yield_se,yield_ci,spectrum=pooled[1]/width,
            spectrum_se=yield_se/width,betas,beta_se,beta_ci,valid_replicates=count(valid),
            invalid_fraction=count(!,valid)/length(valid),replicates)
end

"""
    pbasex_shots(checkpoint; basis, windows, center=nothing, rotation=0, options...)

Estimate integrated electron yields and anisotropies from a schema-1 corrected
event checkpoint with a whole-shot bootstrap. Reads events and SVD blocks from disk.

# Arguments
- `basis` is an existing HDF5 basis. Its geometry selects the folding convention.
  Supply `geometry` to require a matching basis.
- `windows` gives energy intervals with E = alpha × squared detector radius.
- `center` optionally overrides the stored center in unrotated one-based bins.
- `rotation` rotates the histogram and center counterclockwise before folding,
  in degrees divisible by 90. The earlier `quarter_turns` keyword is also supported.
- `shot_weights` gives fixed multipliers in roster order. `run_ids` identifies
  consecutive runs. Set `block_length` for within-run moving-block resampling.
- `replicates`, `confidence`, and `rng` control bootstrap inference.
- `window_batch`, `basis_chunk_bytes`, and `event_chunk` bound read/work batches.
- Optional `replicate_sink(index, draws)` consumes one window's aligned draws
  at a time, for example to write them to disk.

# Returns
Shot count, Legendre orders, bootstrap settings and a summary per window.
Set `keep_replicates=true` to retain aligned yield and beta draws for contrasts.
"""
function pbasex_shots(checkpoint; basis, windows, geometry=nothing, center=nothing, rotation=nothing, quarter_turns::Union{Int,Nothing}=nothing,
        alpha=1.0, reg=0.0, quadrant_filter=[1,1,1,1], shot_weights=nothing,
        run_ids=nothing, block_length::Int=1, replicates::Int=1000, confidence=0.95,
        rng::AbstractRNG=default_rng(), keep_replicates=false, replicate_sink=nothing, window_batch::Int=1,
        basis_chunk_bytes::Int=64*1024^2, event_chunk::Int=65536)
    isfinite(alpha) && alpha > 0 && isfinite(reg) && reg >= 0 || throw(ArgumentError("Invalid alpha or reg"))
    replicates >= 2 && 0 < confidence < 1 && min(block_length,window_batch,basis_chunk_bytes,event_chunk)>0 ||
        throw(ArgumentError("Invalid bootstrap or chunk settings"))
    filter = collect(Int,quadrant_filter)
    length(filter)==4 && all(in((0,1)),filter) && sum(filter)>0 || throw(ArgumentError("Invalid quadrant filter"))
    limits = collect(windows)
    isempty(limits) && throw(ArgumentError("Supply at least one energy interval"))
    isnothing(rotation) || isnothing(quarter_turns) ||
        throw(ArgumentError("Specify only one of rotation and quarter_turns"))
    turns = isnothing(quarter_turns) ? rotation_turns(something(rotation, 0)) : mod(quarter_turns, 4)
    seed = rand(rng,UInt64)
    h5open(checkpoint,"r") do file
        layout = checkpoint_layout(file,turns,center)
        ids = read(file,"shots/id")
        ndims(ids)==1 && eltype(ids)<:Integer && length(ids)>=2 && allunique(ids) ||
            throw(ArgumentError("Expected at least two distinct shots in acquisition order"))
        index = Dict(id=>i for (i,id) in enumerate(ids))
        weights = isnothing(shot_weights) ? ones(length(ids)) : collect(Float64,shot_weights)
        length(weights)==length(ids) && all(w -> isfinite(w) && w>=0,weights) ||
            throw(ArgumentError("Invalid fixed shot weights"))
        starts = [1]
        if !isnothing(run_ids)
            length(run_ids)==length(ids) || throw(ArgumentError("Missing run IDs"))
            append!(starts,(i for i in 2:length(ids) if run_ids[i]!=run_ids[i-1]))
        end
        push!(starts,length(ids)+1)
        runs = [starts[i]:starts[i+1]-1 for i in 1:length(starts)-1]
        all(r -> length(r)>=block_length,runs) || throw(ArgumentError("Block longer than run"))
        any(r -> length(r)>block_length,runs) ||
            throw(ArgumentError("Every run has only one possible block; the bootstrap cannot estimate uncertainty"))
        h5open(basis,"r") do f
            basis_geometry = haskey(f,"geometry") ? Symbol(read(f,"geometry")) : :quadrant
            isnothing(geometry) || geometry === basis_geometry ||
                throw(ArgumentError("Requested geometry $geometry differs from basis geometry $basis_geometry"))
            matlab = !haskey(f,"U") && haskey(f,"Up")
            matlab && basis_geometry !== :quadrant &&
                throw(ArgumentError("MATLAB basis files require quadrant geometry"))
            U = f[matlab ? "Up" : "U"]
            S = vec(Float64.(read(f,"S")))
            V = Float64.(read(f,"V"))
            matlab && (V=permutedims(V))
            x = vec(Float64.(read(f,"x")))
            axial_coordinates(x,basis_geometry)
            orders = vec(Int.(read(f,"l")))
            validate_angular_orders(orders,basis_geometry)
            frk = Float64.(read(f,"frk"))
            scale = alpha/layout.resolution^2
            energy = scale .* x.^2
            all(>(0),diff(energy)) || throw(ArgumentError("Expected increasing energy grid"))
            nk = length(f["k"])
            nl = length(orders)
            pixels = prod(detector_dimensions(length(x),basis_geometry))
            size(frk) == (nk,length(x)) && 0 < length(S) <= min(pixels,nk*nl) &&
                size(U) == (pixels,length(S)) && size(V) == (nk*nl,length(S)) ||
                throw(ArgumentError("Inconsistent basis dimensions"))
            radial = frk .* transpose(x ./ (2scale))
            lookup,fraction,coverage = fold_lookup(layout.shape,layout.center,length(x),basis_geometry,filter)
            lookup,fraction = rotr90(lookup,turns),rotr90(fraction,turns)
            saved_mask = haskey(f,"weights") ? read(f,"weights") : nothing
            isequal(coverage,saved_mask) || throw(ArgumentError("Basis coverage differs from shot geometry"))
            h = inverse_filter(S,reg,pixels,nk*nl)
            # Julia pixel rows are contiguous in native and CPBASEX files.
            columns = max(1,basis_chunk_bytes ÷ (8*pixels))
            summaries = NamedTuple[]
            for first_window in 1:window_batch:length(limits)
                batch = limits[first_window:min(end,first_window+window_batch-1)]
                coordinates = zeros(length(S),nl*length(batch))
                for (i,w) in enumerate(batch)
                    _, integrals = window_integrals(energy,w)
                    q = radial*integrals
                    for j in 1:nl
                        rows = (j-1)*nk+1:j*nk
                        coordinates[:,(i-1)*nl+j] .= h .* (view(V,rows,:)'*q)
                    end
                end
                responses = zeros(pixels,size(coordinates,2))
                for start in 1:columns:length(S)
                    cols = start:min(length(S),start+columns-1)
                    mul!(responses,U[:,cols],view(coordinates,cols,:),1.0,1.0)
                end
                totals = shot_totals(file,layout,ids,index,weights,lookup,fraction,responses,event_chunk)
                samples = resample_totals(totals,runs,block_length,replicates,seed)
                pooled = vec(sum(totals;dims=2))
                for (i,w) in enumerate(batch)
                    rows = (i-1)*nl+1:i*nl
                    summary = window_summary(Tuple(w),view(pooled,rows),view(samples,rows,:),
                                             confidence,keep_replicates || !isnothing(replicate_sink))
                    if !isnothing(replicate_sink)
                        replicate_sink(first_window+i-1,summary.replicates)
                        keep_replicates || (summary=merge(summary,(replicates=nothing,)))
                    end
                    push!(summaries,summary)
                end
            end
            return (; windows=summaries,legendre_orders=orders,shots=length(ids),geometry=basis_geometry,
                    confidence,bootstrap_replicates=replicates,block_length)
        end
    end
end

end
