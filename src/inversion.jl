module PBASEX

using LinearAlgebra

using ..PBASEXTypes: BASISData, PBASEXResult, basis_grid
using ..ImageProcessing: foldQuadrant, resizeFolded, unfoldQuadrant

export pbasex

"""
    pbasex(image; basis=nothing, options...) -> PBASEXResult

Invert a square floating-point detector image with reflection symmetry about
both image axes. The angular orders must start with zero and be even.
Use center=(column, row) and quadrant_filter to select the input pixels.
Selected quadrants are normalized to the same scale as the four-quadrant sum.

Pass an existing Julia or MATLAB CPBASEX HDF5 basis file to use its grid and
precomputed matrices. An existing directory stores bases by definition.
Select `radial_basis=:cubic_bspline` for a compact radial basis when building
a new cache. The default is Gaussian. With basis=nothing, use the package's
data directory. Pass make_images=false
to calculate only spectra, anisotropies, and coefficients.
"""
function pbasex(image::AbstractMatrix{T};
                basis::Union{AbstractString,Nothing}=nothing,
                options...) where {T<:AbstractFloat}
    data = pre_process_image(image, basis, options)
    basis_data = BASISData{T}(data.path; data.basis_options..., data.make_images)
    x_basis, k_basis, l_basis = basis_data.x, basis_data.k, basis_data.l
    validate_angular_orders(l_basis)
    nx, nk, nl = length(x_basis), length(k_basis), length(l_basis)

    U, S, V = basis_data.U, basis_data.S, basis_data.V
    inverse_factors = if iszero(data.reg)
        tolerance = min(size(U, 1), size(V, 1)) * eps(T) * maximum(S)
        map(s -> s > tolerance ? inv(s) : zero(T), S)
    else
        S ./ (S.^2 .+ data.reg)
    end
    coefficients = V * (inverse_factors .* (U' * vec(data.image)))
    ckl = reshape(coefficients, nk, nl)

    radii = T.(x_basis)
    E = data.alpha .* radii.^2
    radial_profiles = basis_data.frk' * ckl
    scaled = (radii ./ (2 * data.alpha)) .* radial_profiles
    IE = scaled[:, 1]
    betas = scaled[:, 2:end] ./ IE
    betas[iszero.(IE), :] .= T(NaN)

    reconstructed_image = inverted_image = residuals = nothing
    if data.make_images
        projection = U * (S .* (V' * coefficients))
        if !isnothing(basis_data.weights)
            projection ./= replace(vec(basis_data.weights), zero(T) => T(NaN))
        end
        reconstructed_image = unfoldQuadrant(reshape(projection, nx, nx))
        slice = (basis_data.Ginv::Matrix{T}) * coefficients
        inverted_image = unfoldQuadrant(reshape(slice, nx, nx))
        if size(reconstructed_image) == size(image) && data.center == (nx, nx)
            residuals = image .- reconstructed_image
        else
            @warn "Residuals require the input and reconstructed detector grids to have the same size and center."
        end
    end
    return PBASEXResult{T}(E, IE, betas, reconstructed_image, inverted_image, residuals, ckl)
end

function validate_angular_orders(l)
    !isempty(l) && first(l) == 0 && allunique(l) && all(order -> order >= 0 && iseven(order), l) ||
        throw(ArgumentError("Image inversion requires distinct, nonnegative, even Legendre orders starting with zero"))
end

function pre_process_image(image::AbstractMatrix{T}, basis_file, options) where T
    rows, cols = size(image)
    @assert rows == cols "Input image must be square, got dimensions ($rows, $cols)"
    basis_keys = (:x, :k, :l, :radial_basis, :params, :rBF, :zIP, :cache_key)
    processing_keys = (:weights, :make_images, :reg, :center, :quadrant_filter, :alpha)
    unknown = setdiff(keys(options), (basis_keys..., processing_keys...))
    isempty(unknown) || @warn "Ignoring unknown keyword arguments" unknown

    make_images = get(options, :make_images, true)::Bool
    reg = T(get(options, :reg, zero(T)))
    alpha = T(get(options, :alpha, one(T)))
    center = get(options, :center, (div(cols, 2) + 1, div(rows, 2) + 1))::Tuple{Int,Int}
    filter = get(options, :quadrant_filter, [1, 1, 1, 1])::Vector{Int}
    weights = get(options, :weights, nothing)::Union{AbstractMatrix{T},Nothing}
    @assert isfinite(reg) && reg >= 0 "Regularization must be finite and nonnegative"
    @assert isfinite(alpha) && alpha > 0 "Energy scaling must be finite and positive"
    if !isnothing(weights)
        @assert size(weights) == size(image) "Input weights must match the image size"
        @assert all(isfinite, weights) "Input weights must be finite"
    end
    haskey(options, :l) && validate_angular_orders(options[:l])

    folded = foldQuadrant(image, center...; quadrant_filter=filter)
    path = isnothing(basis_file) ? mkpath(joinpath(dirname(@__DIR__), "data")) : basis_file
    x = if haskey(options, :x)
        collect(Int, options[:x])
    elseif isfile(path)
        basis_grid(path)
    else
        collect(0:511)
    end
    nx = length(x)
    # Preserve MATLAB's scale for a four-quadrant sum when fewer are selected.
    normalization = T(4) / sum(filter)
    fitted_image = resizeFolded(folded, nx) .* normalization
    fitted_weights = if isnothing(weights)
        nothing
    else
        folded_weights = foldQuadrant(weights, center...; quadrant_filter=filter)
        resizeFolded(folded_weights, nx) .* normalization
    end
    if !isnothing(fitted_weights)
        fitted_image .*= fitted_weights
    end

    basis_options = Dict{Symbol,Any}(key => options[key] for key in basis_keys if haskey(options, key))
    basis_options[:x] = x
    basis_options[:weights] = fitted_weights
    return (; image=fitted_image, path, basis_options, make_images, reg, alpha, center)
end

end # module PBASEX
