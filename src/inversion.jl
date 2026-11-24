module PBASEX

using LinearAlgebra

using ..PBASEXTypes: BASISData, PBASEXResult, basis_grid
using ..ImageProcessing: foldQuadrant, resizeFolded, unfoldQuadrant

export pbasex

"""
    pbasex(image; basis=nothing, options...) -> PBASEXResult

Invert a square floating-point detector image with reflection symmetry about
both image axes. Align the polarization axis vertically along the detector's
`y` coordinate. The Legendre orders must start with zero and be even.
Use `center=(column, row)` and `quadrant_filter` to select the input pixels.
Selected quadrants are normalized to the same scale as the four-quadrant sum.

Pass an existing Julia or MATLAB CPBASEX HDF5 basis file to use its grid and
precomputed matrices. An existing directory stores bases by definition.
Select `radial_basis=:cubic_bspline` for a compact radial basis when building
a new cache. The default is Gaussian. With `basis=nothing`, use the package's
data directory. Pass `make_images=false` to calculate only spectra, anisotropies,
and coefficients.
"""
function pbasex(image::AbstractMatrix{T};
                basis::Union{AbstractString,Nothing}=nothing,
                options...) where {T<:AbstractFloat}
    fit_input = pre_process_image(image, basis, options)
    basis_data = BASISData{T}(fit_input.basis_path; fit_input.basis_options..., fit_input.make_images)
    detector_coordinates, radial_centers, legendre_orders = basis_data.x, basis_data.k, basis_data.l
    validate_angular_orders(legendre_orders)
    detector_side = length(detector_coordinates)
    radial_count = length(radial_centers)
    angular_count = length(legendre_orders)

    U, S, V = basis_data.U, basis_data.S, basis_data.V
    singular_value_filter = if iszero(fit_input.regularization)
        singular_value_cutoff = min(size(U, 1), size(V, 1)) * eps(T) * maximum(S)
        map(singular_value -> singular_value > singular_value_cutoff ? inv(singular_value) : zero(T), S)
    else
        S ./ (S.^2 .+ fit_input.regularization)
    end
    coefficient_vector = V * (singular_value_filter .* (U' * vec(fit_input.fit_pixels)))
    polar_coefficients = reshape(coefficient_vector, radial_count, angular_count)

    radii_pixels = T.(detector_coordinates)
    energy = fit_input.energy_scale .* radii_pixels.^2
    radial_profiles = basis_data.frk' * polar_coefficients
    # Convert each velocity-space profile to a density per unit energy.
    # The radial measure r² dr and E = alpha * r² give r / (2 * alpha).
    energy_profiles = (radii_pixels ./ (2 * fit_input.energy_scale)) .* radial_profiles
    energy_spectrum = energy_profiles[:, 1]
    # Normalize each anisotropic component by the isotropic spectrum.
    anisotropy_parameters = energy_profiles[:, 2:end] ./ energy_spectrum
    anisotropy_parameters[iszero.(energy_spectrum), :] .= T(NaN)

    reconstructed_image = inverted_image = residuals = nothing
    if fit_input.make_images
        projected_pixels = U * (S .* (V' * coefficient_vector))
        if !isnothing(basis_data.weights)
            projected_pixels ./= replace(vec(basis_data.weights), zero(T) => T(NaN))
        end
        reconstructed_image = unfoldQuadrant(reshape(projected_pixels, detector_side, detector_side))
        slice_pixels = (basis_data.Ginv::Matrix{T}) * coefficient_vector
        inverted_image = unfoldQuadrant(reshape(slice_pixels, detector_side, detector_side))
        if size(reconstructed_image) == size(image) && fit_input.center == (detector_side, detector_side)
            residuals = image .- reconstructed_image
        else
            @warn "Residuals require the input and reconstructed detector grids to have the same size and center."
        end
    end
    return PBASEXResult{T}(energy, energy_spectrum, anisotropy_parameters,
                           reconstructed_image, inverted_image, residuals, polar_coefficients)
end

function validate_angular_orders(legendre_orders)
    !isempty(legendre_orders) && first(legendre_orders) == 0 && allunique(legendre_orders) && all(order -> order >= 0 && iseven(order), legendre_orders) ||
        throw(ArgumentError("Image inversion requires distinct, nonnegative, even Legendre orders starting with zero"))
end

function pre_process_image(image::AbstractMatrix{T}, basis_file, options) where T
    image_rows, image_columns = size(image)
    @assert image_rows == image_columns "Input image must be square, got dimensions ($image_rows, $image_columns)"
    basis_keys = (:x, :k, :l, :radial_basis, :params, :rBF, :zIP, :cache_key)
    processing_keys = (:weights, :make_images, :reg, :center, :quadrant_filter, :alpha)
    unknown_options = setdiff(keys(options), (basis_keys..., processing_keys...))
    isempty(unknown_options) || @warn "Ignoring unknown keyword arguments" unknown=unknown_options

    make_images = get(options, :make_images, true)::Bool
    regularization = T(get(options, :reg, zero(T)))
    energy_scale = T(get(options, :alpha, one(T)))
    center = get(options, :center, (div(image_columns, 2) + 1, div(image_rows, 2) + 1))::Tuple{Int,Int}
    quadrant_filter = get(options, :quadrant_filter, [1, 1, 1, 1])::Vector{Int}
    weights = get(options, :weights, nothing)::Union{AbstractMatrix{T},Nothing}
    @assert isfinite(regularization) && regularization >= 0 "Regularization must be finite and nonnegative"
    @assert isfinite(energy_scale) && energy_scale > 0 "Energy scaling must be finite and positive"
    if !isnothing(weights)
        @assert size(weights) == size(image) "Input weights must match the image size"
        @assert all(isfinite, weights) "Input weights must be finite"
    end
    haskey(options, :l) && validate_angular_orders(options[:l])

    folded_projection = foldQuadrant(image, center...; quadrant_filter=quadrant_filter)
    basis_path = isnothing(basis_file) ? mkpath(joinpath(dirname(@__DIR__), "data")) : basis_file
    detector_coordinates = if haskey(options, :x)
        collect(Int, options[:x])
    elseif isfile(basis_path)
        basis_grid(basis_path)
    else
        collect(0:511)
    end
    detector_side = length(detector_coordinates)
    # Preserve MATLAB's scale for a four-quadrant sum when fewer are selected.
    quadrant_normalization = T(4) / sum(quadrant_filter)
    fit_pixels = resizeFolded(folded_projection, detector_side) .* quadrant_normalization
    fit_weights = if isnothing(weights)
        nothing
    else
        folded_pixel_weights = foldQuadrant(weights, center...; quadrant_filter=quadrant_filter)
        resizeFolded(folded_pixel_weights, detector_side) .* quadrant_normalization
    end
    # Weight the measured projection consistently with the forward operator.
    if !isnothing(fit_weights)
        fit_pixels .*= fit_weights
    end

    basis_options = Dict{Symbol,Any}(key => options[key] for key in basis_keys if haskey(options, key))
    basis_options[:x] = detector_coordinates
    basis_options[:weights] = fit_weights
    return (; fit_pixels, basis_path, basis_options, make_images, regularization, energy_scale, center)
end

end # module PBASEX
