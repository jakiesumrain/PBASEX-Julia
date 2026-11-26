module PBASEXTypes

using HDF5
using FileWatching.Pidfile: trymkpidlock

using ..BasisFunctions: build_forward_basis, cubic_bspline, cubic_bspline_limit
using ..InverseBasisFunctions: build_central_slice_basis
using ..MathUtils: meshgrid, detector_dimensions, axial_coordinates

export PBASEXResult, BASISData

"""
    PBASEXResult{T}

Results for one VMI image. `E` contains energy values on the requested scale,
and `IE` contains the one-dimensional energy spectrum. Columns of
`betas` follow the nonzero `legendre_orders`. The coefficient
matrix has one row per radial center and one column per Legendre order.
`inverted_image` is the central velocity slice. `reconstructed_image` is
the fitted detector projection. Image outputs and residuals are optional.
"""
struct PBASEXResult{T<:AbstractFloat}
    E::Vector{T}
    IE::Vector{T}
    betas::Matrix{T}
    reconstructed_image::Union{Matrix{T},Nothing}
    inverted_image::Union{Matrix{T},Nothing}
    residuals::Union{Matrix{T},Nothing}
    coefficients::Matrix{T}
    legendre_orders::Union{Vector{Int},Nothing}
end

PBASEXResult{T}(E, IE, betas, reconstructed_image, inverted_image, residuals, coefficients) where {T<:AbstractFloat} =
    PBASEXResult{T}(E, IE, betas, reconstructed_image, inverted_image, residuals, coefficients, nothing)

PBASEXResult(E::Vector{T}, IE, betas, reconstructed_image, inverted_image, residuals, coefficients) where {T<:AbstractFloat} =
    PBASEXResult{T}(E, IE, betas, reconstructed_image, inverted_image, residuals, coefficients)

"""
    BASISData{T}(; x, k, l, radial_basis=:gaussian, geometry=:quadrant, params, weights, rBF, zIP, make_images=true)
    BASISData{T}(path; cache_key=nothing, kwargs...)
    BASISData([path]; kwargs...)

Compute a basis in memory or supply a file path to load or save it.
If the file exists, any specified basis settings must match its metadata.
The fields `x`, `k`, and `l` give the detector coordinates, radial basis
centers, and Legendre orders. The columns of `U` and `V` are the left and
right singular vectors of the forward operator. `S` holds the singular
values. The radial table `frk` has one row per radial center and one column
per sampled radius. `Ginv` maps coefficients to a central slice and is
omitted if `make_images=false`. Existing MATLAB CPBASEX files are supported.
Use `geometry=:half_plane` to retain both axial hemispheres and fit odd
Legendre orders.

Select `radial_basis=:cubic_bspline` for a compact cubic B-spline. Its first
parameter is the knot spacing in pixels. For a custom `rBF` or `zIP`,
`cache_key` records an optional identity string. Supply the same string
if you pass the original functions again when opening the file.
"""
struct BASISData{T<:AbstractFloat}
    x::Vector{Int}
    k::Vector{T}
    l::Vector{Int}
    weights::Union{Matrix{T},Nothing}
    params::Vector{T}
    rBF::Function
    zIP::Function
    U::Matrix{T}
    S::Vector{T}
    V::Matrix{T}
    frk::Matrix{T}
    Ginv::Union{Matrix{T},Nothing}
    geometry::Symbol
end

BASISData{T}(x, k, l, weights, params, rBF, zIP, U, S, V, frk, Ginv) where {T<:AbstractFloat} =
    BASISData{T}(x, k, l, weights, params, rBF, zIP, U, S, V, frk, Ginv, :quadrant)

BASISData(x, k::Vector{T}, l, weights, params, rBF, zIP, U, S, V, frk, Ginv) where {T<:AbstractFloat} =
    BASISData{T}(x, k, l, weights, params, rBF, zIP, U, S, V, frk, Ginv)

gaussian(radius, radial_center, basis_parameters) = exp(-(radius - radial_center)^2 / (2 * basis_parameters[1]^2))
function gaussian_limit(detector_radius, radial_center::T, basis_parameters) where T
    # Keep the CPBASEX cutoff at exp(-10) of the Gaussian peak.
    return sqrt(max(zero(T), (sqrt(T(20)) * basis_parameters[1] + radial_center)^2 - detector_radius^2))
end

function radial_functions(family)
    family === :gaussian && return gaussian, gaussian_limit
    family === :cubic_bspline && return cubic_bspline, cubic_bspline_limit
    throw(ArgumentError("radial_basis must be :gaussian or :cubic_bspline"))
end

function radial_family(data)
    return data.rBF === cubic_bspline && data.zIP === cubic_bspline_limit ? :cubic_bspline : :gaussian
end

function basis_definition(::Type{T}, options) where T
    detector_coordinates = collect(Int, get(options, :x, 0:511))
    isempty(detector_coordinates) && throw(ArgumentError("x must contain at least one grid point"))
    geometry = get(options, :geometry, :quadrant)
    grid_size = detector_dimensions(length(detector_coordinates), geometry)
    axial_coordinates(detector_coordinates, geometry)
    radial_centers = collect(T, get(options, :k, detector_coordinates[2:4:end] .+ 0.5))
    isempty(radial_centers) && throw(ArgumentError("k must contain at least one radial center"))
    default_orders = geometry === :half_plane ? (0, 1, 2, 3, 4) : (0, 2, 4)
    legendre_orders = collect(Int, get(options, :l, default_orders))
    !isempty(legendre_orders) && all(>=(0), legendre_orders) || throw(ArgumentError("l must contain nonnegative orders"))
    basis_family = get(options, :radial_basis, :gaussian)
    default_radial_profile, default_line_of_sight_limit = radial_functions(basis_family)
    default_width = basis_family === :gaussian ?
        0.7 * length(detector_coordinates) / length(radial_centers) :
        length(detector_coordinates) / length(radial_centers)
    basis_parameters = collect(T, get(options, :params, [default_width]))
    pixel_weights = get(options, :weights, nothing)
    if !isnothing(pixel_weights)
        size(pixel_weights) == grid_size ||
            throw(ArgumentError("Basis weights must match the detector grid $grid_size"))
        pixel_weights = Matrix{T}(pixel_weights)
    end
    radial_profile = get(options, :rBF, default_radial_profile)::Function
    line_of_sight_limit = get(options, :zIP, default_line_of_sight_limit)::Function
    if haskey(options, :radial_basis) && (radial_profile !== default_radial_profile || line_of_sight_limit !== default_line_of_sight_limit)
        throw(ArgumentError("Choose a built-in radial_basis or supply custom rBF and zIP"))
    end
    all(isfinite, radial_centers) && all(isfinite, basis_parameters) || throw(ArgumentError("k and params must be finite"))
    if radial_profile === gaussian || line_of_sight_limit === gaussian_limit
        !isempty(basis_parameters) && first(basis_parameters) > 0 ||
            throw(ArgumentError("The Gaussian basis requires a positive width"))
    elseif basis_family === :cubic_bspline
        !isempty(basis_parameters) && first(basis_parameters) > 0 ||
            throw(ArgumentError("The cubic B-spline requires a positive knot spacing"))
    end
    # Use the public keyword names so the definition can be passed to the constructor.
    return (; x=detector_coordinates, k=radial_centers, l=legendre_orders, geometry,
             params=basis_parameters, weights=pixel_weights, rBF=radial_profile, zIP=line_of_sight_limit)
end

function BASISData{T}(; make_images::Bool=true, kwargs...) where {T<:AbstractFloat}
    definition = basis_definition(T, kwargs)
    if all(key -> haskey(kwargs, key), (:U, :S, :V))
        U, S, V = kwargs[:U], kwargs[:S], kwargs[:V]
    else
        U, S, V = build_forward_basis(definition.x, definition.k, definition.l, definition.rBF,
                                      definition.zIP, definition.params, definition.weights; definition.geometry)
    end
    radial_table = get(kwargs, :frk) do
        radius_grid, radial_center_grid = meshgrid(T.(definition.x), definition.k)
        definition.rBF.(radius_grid, radial_center_grid, Ref(definition.params))
    end
    slice_basis = if make_images
        get(kwargs, :Ginv) do
            build_central_slice_basis(definition.x, definition.k, definition.l, definition.rBF,
                                      definition.params; definition.geometry)
        end
    else
        nothing
    end
    return BASISData{T}(definition.x, definition.k, definition.l, definition.weights,
                        definition.params, definition.rBF, definition.zIP, U, S, V,
                        radial_table, slice_basis, definition.geometry)
end

BASISData(; kwargs...) = BASISData{Float64}(; kwargs...)

function BASISData{T}(path::AbstractString; make_images::Bool=true,
                      cache_key::Union{AbstractString,Nothing}=nothing, kwargs...) where {T<:AbstractFloat}
    isdir(path) && throw(ArgumentError("Choose a basis filename inside $path"))
    isfile(path) && return load_basis(T, path, make_images, cache_key, kwargs)
    isdir(dirname(abspath(path))) || throw(ArgumentError("Parent directory of basis file does not exist: $path"))

    data = BASISData{T}(; make_images, kwargs...)
    save_basis_data(data, path, cache_key) && return data
    isfile(path) || throw(ArgumentError("Another process is publishing the basis at $path. Retry once the file is available."))
    return load_basis(T, path, make_images, cache_key, kwargs)
end

BASISData(path::AbstractString; kwargs...) = BASISData{Float64}(path; kwargs...)

function basis_grid(path)
    return h5open(path, "r") do file
        return vec(Int.(read(file, "x")))
    end
end

function unavailable_basis_function(args...)
    throw(ArgumentError("This file does not record the radial basis function. Supply the original function to regenerate a slice."))
end

function load_basis(::Type{T}, path, make_images, cache_key, options) where T
    return h5open(path, "r") do basis_file
        is_matlab_basis = haskey(basis_file, "Up") && !haskey(basis_file, "U")
        detector_coordinates = vec(Int.(read(basis_file, "x")))
        geometry = haskey(basis_file, "geometry") ? Symbol(read(basis_file, "geometry")) : :quadrant
        grid_size = detector_dimensions(length(detector_coordinates), geometry)
        axial_coordinates(detector_coordinates, geometry)
        is_matlab_basis && geometry !== :quadrant &&
            throw(ArgumentError("MATLAB basis files require quadrant geometry"))
        get(options, :geometry, geometry) === geometry ||
            throw(ArgumentError("Requested geometry does not match the basis in $path"))
        radial_centers = vec(T.(read(basis_file, "k")))
        legendre_orders = vec(Int.(read(basis_file, "l")))
        pixel_weights = haskey(basis_file, "weights") ? convert(Matrix{T}, read(basis_file, "weights")) : nothing
        has_basis_parameters = haskey(basis_file, "params")
        basis_parameters = has_basis_parameters ? vec(T.(read(basis_file, "params"))) : T[]
        basis_family = haskey(basis_file, "radial_basis") ? Symbol(read(basis_file, "radial_basis")) : :gaussian
        default_radial_profile, default_line_of_sight_limit = radial_functions(basis_family)
        custom_functions = haskey(basis_file, "custom_functions") ?
            Tuple(Bool.(read(basis_file, "custom_functions"))) :
            (!has_basis_parameters, !has_basis_parameters)
        saved_key = haskey(basis_file, "cache_key") ? read(basis_file, "cache_key") : nothing

        if haskey(options, :radial_basis) && (options[:radial_basis] !== basis_family || any(custom_functions))
            throw(ArgumentError("Requested radial_basis does not match the basis in $path"))
        end
        for (name, value) in pairs((; x=detector_coordinates, k=radial_centers, l=legendre_orders))
            if haskey(options, name) && !isequal(collect(eltype(value), options[name]), value)
                throw(ArgumentError("Requested $name does not match the basis in $path"))
            end
        end
        if haskey(options, :weights)
            requested_weights = options[:weights]
            requested_weights = isnothing(requested_weights) ? nothing : Matrix{T}(requested_weights)
            isequal(requested_weights, pixel_weights) ||
                throw(ArgumentError("Requested weights do not match the basis in $path"))
        end
        if haskey(options, :params) && (!has_basis_parameters || !isequal(collect(T, options[:params]), basis_parameters))
            throw(ArgumentError("Requested params could not be verified against the basis in $path"))
        end
        if !isnothing(cache_key) && cache_key != saved_key
            throw(ArgumentError("Requested cache_key does not match the basis in $path"))
        end
        radial_profile = get(options, :rBF, custom_functions[1] ? unavailable_basis_function : default_radial_profile)
        line_of_sight_limit = get(options, :zIP, custom_functions[2] ? unavailable_basis_function : default_line_of_sight_limit)
        for (name, default, was_custom) in zip((:rBF, :zIP), (default_radial_profile, default_line_of_sight_limit), custom_functions)
            if haskey(options, name)
                if (options[name] !== default) != was_custom || (was_custom && isnothing(cache_key))
                    throw(ArgumentError("Cannot verify $name against the basis in $path. Use its matching cache_key or load the saved functions' matrices by path."))
                end
            end
        end
        if !is_matlab_basis && !has_basis_parameters
            @warn "Legacy Julia basis has no weight metadata. Loading it as unweighted." path
        end

        # In a MATLAB CPBASEX cache, Up already has the shape of Julia's U.
        # The stored V and central-slice matrix need transposing.
        U = convert(Matrix{T}, read(basis_file, is_matlab_basis ? "Up" : "U"))
        S = vec(T.(read(basis_file, "S")))
        V = convert(Matrix{T}, read(basis_file, "V"))
        radial_table = convert(Matrix{T}, read(basis_file, "frk"))
        is_matlab_basis && (V = permutedims(V))
        pixel_count = prod(grid_size)
        basis_count = length(radial_centers) * length(legendre_orders)
        singular_value_count = length(S)
        size(U) == (pixel_count, singular_value_count) && size(V) == (basis_count, singular_value_count) &&
            size(radial_table) == (length(radial_centers), length(detector_coordinates)) && 0 < singular_value_count <= min(pixel_count, basis_count) &&
            (isnothing(pixel_weights) || size(pixel_weights) == grid_size) ||
            throw(ArgumentError("Basis dimensions are inconsistent in $path"))
        slice_basis = nothing
        if make_images
            if haskey(basis_file, "Ginv")
                slice_basis = convert(Matrix{T}, read(basis_file, "Ginv"))
                is_matlab_basis && (slice_basis = permutedims(slice_basis))
            else
                radial_profile === unavailable_basis_function &&
                    throw(ArgumentError("The basis in $path has no central-slice matrix or radial function"))
                slice_basis = build_central_slice_basis(detector_coordinates, radial_centers,
                                                        legendre_orders, radial_profile, basis_parameters; geometry)
            end
            size(slice_basis) == (pixel_count, basis_count) ||
                throw(ArgumentError("Central-slice matrix has inconsistent dimensions in $path"))
        end
        return BASISData{T}(detector_coordinates, radial_centers, legendre_orders, pixel_weights,
                            basis_parameters, radial_profile, line_of_sight_limit, U, S, V,
                            radial_table, slice_basis, geometry)
    end
end

function save_basis_data(data::BASISData, path::AbstractString, cache_key)
    # Write in the destination directory so a reader never opens a partial file.
    return mktemp(dirname(abspath(path))) do temp, io
        close(io)
        h5open(temp, "w") do f
            for name in (:x, :k, :l, :params, :weights, :U, :S, :V, :frk, :Ginv)
                value = getproperty(data, name)
                isnothing(value) || write(f, String(name), value)
            end
            family = radial_family(data)
            default_rBF, default_zIP = radial_functions(family)
            write(f, "custom_functions", UInt8[data.rBF !== default_rBF, data.zIP !== default_zIP])
            write(f, "radial_basis", String(family))
            write(f, "geometry", String(data.geometry))
            isnothing(cache_key) || write(f, "cache_key", String(cache_key))
        end
        # Serialize publication when different processes choose the same filename.
        trymkpidlock(abspath(path) * ".lock") do
            isfile(path) && return false
            mv(temp, path)
            return true
        end
    end
end

end # module PBASEXTypes
