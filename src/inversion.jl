module Inversion

using LinearAlgebra

using ..PBASEXTypes: BASISData, PBASEXResult, basis_grid
using ..ImageProcessing: fold_quadrants, resize_folded_quadrant, unfold_quadrant,
                         fold_half_plane, unfold_half_plane, quadrant_sum
using ..CountingErrors: summarize_counting_errors
using ..MathUtils: detector_dimensions, axial_coordinates

export pbasex, pbasex_counts

"""
    pbasex(image; basis, geometry=:quadrant, options...) -> PBASEXResult

Invert a floating-point detector image with a vertical physical symmetry
axis. The default `geometry=:quadrant` folds about both image axes and fits
even Legendre orders. Set `geometry=:half_plane` to retain both axial
hemispheres and fit odd orders as well. In half-plane mode, positive axial
coordinates point toward the top of the image.

Pass a Julia or MATLAB CPBASEX HDF5 file path as `basis`. An existing
file supplies the grid and matrices. A new path creates a file using
the requested basis settings. Select `radial_basis=:cubic_bspline` to
build a compact radial basis. The default is Gaussian. Pass
`make_images=false` to calculate only spectra, anisotropies, and coefficients.

# Arguments
- `image` holds measured detector intensities.
- `basis` is the path of the HDF5 basis file.
- `geometry` selects quadrant or half-plane fitting and must match the basis file.
- `center=(column, row)` and `quadrant_filter` select the detector pixels.

# Returns
A `PBASEXResult` with the fitted spectrum, anisotropies and their Legendre
orders, coefficients, and optional reconstructed images.
"""
function pbasex(image::AbstractMatrix{T};
                basis::AbstractString,
                geometry=:quadrant,
                options...) where {T<:AbstractFloat}
    return fit_projection(image, basis, geometry, options, nothing)
end

"""
    pbasex_counts(image; basis, geometry=:quadrant, peak_windows=(), options...)

Fit raw detector counts and propagate independent Poisson counting noise.

# Arguments
- `image` is a nonnegative, integer-valued detector histogram.
- `peak_windows` contains pairs of energy limits in the units of `result.E`.
- `basis`, `geometry` and fit options follow `pbasex`. Use a fixed center and
  the original counts before interpolation, subtraction or intensity scaling.

# Returns
A named tuple with the ordinary fit (`result`), pointwise one-standard-error
bars (`spectrum_se`), and one summary per window (`peaks`). Each summary
contains `window`, `yield`, `yield_se`, `betas`, and `beta_se`. Beta columns
follow `result.legendre_orders[2:end]`. A yield within one standard error
of zero or below gives undefined symmetric beta error bars (`NaN`).
"""
function pbasex_counts(image::AbstractMatrix{<:Real}; basis::AbstractString,
                       geometry=:quadrant, peak_windows=(), options...)
    all(n -> isfinite(n) && n >= 0 && isinteger(n), image) ||
        throw(ArgumentError("Counting errors require finite nonnegative raw integer counts"))
    isnothing(get(options, :weights, nothing)) ||
        throw(ArgumentError("Counting errors require the raw unweighted detector histogram"))
    counts = eltype(image) <: AbstractFloat ? image : Float64.(image)
    return fit_projection(counts, basis, geometry, options, collect(peak_windows))
end

function fit_projection(image::AbstractMatrix{T}, basis, geometry, options, windows) where T
    fit_input = pre_process_image(image, basis, options, geometry; counting=!isnothing(windows))
    basis_data = BASISData{T}(basis; fit_input.basis_options..., fit_input.make_images)
    detector_coordinates, radial_centers, legendre_orders = basis_data.x, basis_data.k, basis_data.l
    validate_angular_orders(legendre_orders, geometry)
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
        grid_size = detector_dimensions(detector_side, geometry)
        unfold = geometry === :half_plane ? unfold_half_plane : mirror_quadrant
        reconstructed_image = unfold(reshape(projected_pixels, grid_size))
        slice_pixels = (basis_data.Ginv::Matrix{T}) * coefficient_vector
        # CPBASEX quadrant slice matrices omit the angular density factor.
        if geometry === :quadrant
            slice_pixels ./= 4 * T(π)
        end
        inverted_image = unfold(reshape(slice_pixels, grid_size))
        if size(reconstructed_image) == size(image) && fit_input.center == (detector_side, detector_side)
            residuals = image .- reconstructed_image
        else
            @warn "Residuals require the input and reconstructed detector grids to have the same size and center."
        end
    end
    result = PBASEXResult{T}(energy, energy_spectrum, anisotropy_parameters,
                           reconstructed_image, inverted_image, residuals, polar_coefficients,
                           legendre_orders)
    isnothing(windows) && return result
    return summarize_counting_errors(result, basis_data, singular_value_filter,
                                     fit_input.count_variance, fit_input.energy_scale, windows)
end

mirror_quadrant(pixels) = unfold_quadrant(pixels; quadrant_count=1)

function validate_angular_orders(legendre_orders, geometry)
    !isempty(legendre_orders) && first(legendre_orders) == 0 &&
        allunique(legendre_orders) && all(>=(0), legendre_orders) ||
        throw(ArgumentError("Image inversion requires distinct, nonnegative Legendre orders starting with zero"))
    geometry === :quadrant && any(isodd, legendre_orders) &&
        throw(ArgumentError("Odd Legendre orders require geometry=:half_plane"))
end

function pre_process_image(image::AbstractMatrix{T}, basis_file, options, geometry=:quadrant;
                           counting=false) where T
    image_rows, image_columns = size(image)
    detector_dimensions(1, geometry)
    @assert geometry === :half_plane || image_rows == image_columns "Input image must be square, got dimensions ($image_rows, $image_columns)"
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
    haskey(options, :l) && validate_angular_orders(options[:l], geometry)

    detector_coordinates = if haskey(options, :x)
        collect(Int, options[:x])
    elseif isfile(basis_file)
        basis_grid(basis_file)
    else
        collect(0:511)
    end
    detector_side = length(detector_coordinates)
    count_variance = nothing
    if geometry === :half_plane
        detector_side > 0 || throw(ArgumentError("x must contain at least one grid point"))
        axial_coordinates(detector_coordinates, geometry)
        counting && (count_variance = zeros(T, 2 * detector_side - 1, detector_side))
        fit_pixels, fit_weights = fold_half_plane(image, center, detector_side, quadrant_filter, weights;
                                                count_variance)
    else
        folded_projection = fold_quadrants(image, center...; quadrant_filter=quadrant_filter)
        selected_count = sum(quadrant_filter)
        fit_pixels = resize_folded_quadrant(folded_projection, detector_side) ./ T(selected_count)
        if counting
            variance_sum = quadrant_sum(image, center..., quadrant_filter, true)
            count_variance = resize_folded_quadrant(variance_sum, detector_side) ./ T(selected_count^2)
        end
        fit_weights = if isnothing(weights)
            nothing
        else
            folded_pixel_weights = fold_quadrants(weights, center...; quadrant_filter=quadrant_filter)
            # Preserve the forward-operator scale recorded in existing basis files.
            resize_folded_quadrant(folded_pixel_weights, detector_side) .* (T(4) / selected_count)
        end
    end
    # Weight the measured projection consistently with the forward operator.
    if !isnothing(fit_weights)
        fit_pixels .*= fit_weights
    end

    basis_options = Dict{Symbol,Any}(key => options[key] for key in basis_keys if haskey(options, key))
    basis_options[:x] = detector_coordinates
    basis_options[:weights] = fit_weights
    basis_options[:geometry] = geometry
    return (; fit_pixels, basis_options, make_images, regularization, energy_scale, center, count_variance)
end

end # module Inversion
