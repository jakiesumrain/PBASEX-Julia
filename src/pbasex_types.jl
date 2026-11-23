module PBASEXTypes

using HDF5
using Serialization
using SHA

using ..BasisFunctions: findG
using ..InverseBasisFunctions: findGinv
using ..MathUtils: meshgrid

export PBASEXResult, BASISData, generate_hashed_filename

"""
    PBASEXResult{T}

Results for one image. E and IE are the energy scale and one-dimensional
intensity spectrum. Each column of betas corresponds to a nonzero Legendre
order in the basis. Coefficients have one row per radial center and one
column per Legendre order. Image outputs and residuals are optional.
"""
struct PBASEXResult{T<:AbstractFloat}
    E::Vector{T}
    IE::Vector{T}
    betas::Matrix{T}
    reconstructed_image::Union{Matrix{T},Nothing}
    inverted_image::Union{Matrix{T},Nothing}
    residuals::Union{Matrix{T},Nothing}
    coefficients::Matrix{T}
end

"""
    BASISData{T}(; x, k, l, params, weights, rBF, zIP, make_images=true)
    BASISData{T}(path; cache_key=nothing, kwargs...)
    BASISData(path; kwargs...)

Compute a basis, cache it in a directory, or load it from an explicit file.
The columns of U and V are the left and right singular vectors of the
projected basis. The radial table frk has one row per radial center.
Ginv maps coefficients to a central slice and is omitted if
make_images=false. Existing MATLAB CPBASEX basis files are supported.

The directory cache includes the grid, radial centers, angular orders,
parameters and pixel weights. Pass a stable cache_key when caching a
custom rBF or zIP in a directory. A custom basis without a key is computed
in memory in directory mode. An explicit file can still store a custom
basis and be reloaded by path.
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
end

gaussian(r, k, p) = exp(-(r - k)^2 / (2 * p[1]^2))
function gaussian_limit(R, k::T, p) where T
    return sqrt(max(zero(T), (sqrt(T(20)) * p[1] + k)^2 - R^2))
end

function basis_definition(::Type{T}, options) where T
    x = collect(Int, get(options, :x, 0:511))
    isempty(x) && throw(ArgumentError("x must contain at least one grid point"))
    k = collect(T, get(options, :k, x[2:4:end] .+ 0.5))
    isempty(k) && throw(ArgumentError("k must contain at least one radial center"))
    l = collect(Int, get(options, :l, (0, 2, 4)))
    !isempty(l) && all(>=(0), l) || throw(ArgumentError("l must contain nonnegative orders"))
    params = collect(T, get(options, :params, [0.7 * length(x) / length(k)]))
    weights = get(options, :weights, nothing)
    if !isnothing(weights)
        size(weights) == (length(x), length(x)) ||
            throw(ArgumentError("Basis weights must match the square grid"))
        weights = Matrix{T}(weights)
    end
    rBF = get(options, :rBF, gaussian)::Function
    zIP = get(options, :zIP, gaussian_limit)::Function
    all(isfinite, k) && all(isfinite, params) || throw(ArgumentError("k and params must be finite"))
    if rBF === gaussian || zIP === gaussian_limit
        !isempty(params) && first(params) > 0 ||
            throw(ArgumentError("The Gaussian basis requires a positive width"))
    end
    return (; x, k, l, params, weights, rBF, zIP)
end

function BASISData{T}(; make_images::Bool=true, kwargs...) where {T<:AbstractFloat}
    d = basis_definition(T, kwargs)
    if all(key -> haskey(kwargs, key), (:U, :S, :V))
        U, S, V = kwargs[:U], kwargs[:S], kwargs[:V]
    else
        U, S, V = findG(d.x, d.k, d.l, d.rBF, d.zIP, d.params, d.weights)
    end
    frk = get(kwargs, :frk) do
        R, K = meshgrid(T.(d.x), d.k)
        d.rBF.(R, K, Ref(d.params))
    end
    Ginv = if make_images
        get(kwargs, :Ginv) do
            findGinv(d.x, d.k, d.l, d.rBF, d.params)
        end
    else
        nothing
    end
    return BASISData{T}(d.x, d.k, d.l, d.weights, d.params, d.rBF, d.zIP,
                        U, S, V, frk, Ginv)
end

function BASISData{T}(path::AbstractString; make_images::Bool=true,
                      cache_key::Union{AbstractString,Nothing}=nothing, kwargs...) where {T<:AbstractFloat}
    options = (; kwargs...)
    if isdir(path)
        d = basis_definition(T, options)
        custom = (d.rBF !== gaussian, d.zIP !== gaussian_limit)
        if any(custom) && isnothing(cache_key)
            return BASISData{T}(; make_images, d...)
        end
        identity = (d.params, d.weights, custom, isnothing(cache_key) ? nothing : String(cache_key))
        filename = hashed_filename((:pbasex_v2, d.x, d.k, d.l, identity))
        path = joinpath(path, filename)
        options = merge(options, d)
    end
    isfile(path) && return load_basis(T, path, make_images, cache_key, options)
    isdir(dirname(abspath(path))) || throw(ArgumentError("Parent directory of basis file does not exist: $path"))

    data = BASISData{T}(; make_images, options...)
    save_basis_data(data, path, cache_key)
    return data
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
    return h5open(path, "r") do f
        matlab = haskey(f, "Up") && !haskey(f, "U")
        x = vec(Int.(read(f, "x")))
        k = vec(T.(read(f, "k")))
        l = vec(Int.(read(f, "l")))
        weights = haskey(f, "weights") ? convert(Matrix{T}, read(f, "weights")) : nothing
        has_params = haskey(f, "params")
        params = has_params ? vec(T.(read(f, "params"))) : T[]
        custom = haskey(f, "custom_functions") ? Tuple(Bool.(read(f, "custom_functions"))) : (!has_params, !has_params)
        saved_key = haskey(f, "cache_key") ? read(f, "cache_key") : nothing

        for (name, value) in pairs((; x, k, l))
            if haskey(options, name) && !isequal(collect(eltype(value), options[name]), value)
                throw(ArgumentError("Requested $name does not match the basis in $path"))
            end
        end
        if haskey(options, :weights)
            requested_weights = options[:weights]
            requested_weights = isnothing(requested_weights) ? nothing : Matrix{T}(requested_weights)
            isequal(requested_weights, weights) ||
                throw(ArgumentError("Requested weights do not match the basis in $path"))
        end
        if haskey(options, :params) && (!has_params || !isequal(collect(T, options[:params]), params))
            throw(ArgumentError("Requested params could not be verified against the basis in $path"))
        end
        if !isnothing(cache_key) && cache_key != saved_key
            throw(ArgumentError("Requested cache_key does not match the basis in $path"))
        end
        rBF = get(options, :rBF, custom[1] ? unavailable_basis_function : gaussian)
        zIP = get(options, :zIP, custom[2] ? unavailable_basis_function : gaussian_limit)
        for (name, default, was_custom) in zip((:rBF, :zIP), (gaussian, gaussian_limit), custom)
            if haskey(options, name)
                if (options[name] !== default) != was_custom || (was_custom && isnothing(cache_key))
                    throw(ArgumentError("Cannot verify $name against the basis in $path. Use its matching cache_key or load the saved functions' matrices by path."))
                end
            end
        end
        if !matlab && !has_params
            @warn "Legacy Julia basis has no weight metadata. Loading it as unweighted." path
        end

        # CPBASEX's loadG transposes Up, V and Ginv from the MATLAB file.
        # The Julia solver stores U (rather than U') and uses frk[k, x].
        U = convert(Matrix{T}, read(f, matlab ? "Up" : "U"))
        S = vec(T.(read(f, "S")))
        V = convert(Matrix{T}, read(f, "V"))
        frk = convert(Matrix{T}, read(f, "frk"))
        matlab && (V = permutedims(V))
        pixels = length(x)^2
        columns = length(k) * length(l)
        rank_size = length(S)
        size(U) == (pixels, rank_size) && size(V) == (columns, rank_size) &&
            size(frk) == (length(k), length(x)) && 0 < rank_size <= min(pixels, columns) &&
            (isnothing(weights) || size(weights) == (length(x), length(x))) ||
            throw(ArgumentError("Basis dimensions are inconsistent in $path"))
        Ginv = nothing
        if make_images
            if haskey(f, "Ginv")
                Ginv = convert(Matrix{T}, read(f, "Ginv"))
                matlab && (Ginv = permutedims(Ginv))
            else
                rBF === unavailable_basis_function &&
                    throw(ArgumentError("The basis in $path has no central-slice matrix or radial function"))
                Ginv = findGinv(x, k, l, rBF, params)
            end
            size(Ginv) == (pixels, columns) ||
                throw(ArgumentError("Central-slice matrix has inconsistent dimensions in $path"))
        end
        return BASISData{T}(x, k, l, weights, params, rBF, zIP, U, S, V, frk, Ginv)
    end
end

function save_basis_data(data::BASISData, path::AbstractString, cache_key)
    # Write in the destination directory so a reader never opens a partial file.
    mktemp(dirname(abspath(path))) do temp, io
        close(io)
        h5open(temp, "w") do f
            for name in (:x, :k, :l, :params, :weights, :U, :S, :V, :frk, :Ginv)
                value = getproperty(data, name)
                isnothing(value) || write(f, String(name), value)
            end
            write(f, "custom_functions", UInt8[data.rBF !== gaussian, data.zIP !== gaussian_limit])
            isnothing(cache_key) || write(f, "cache_key", String(cache_key))
        end
        try
            mv(temp, path)
        catch err
            isfile(path) || rethrow(err)
            # A competing writer already published a complete cache.
        end
    end
end

"""
    generate_hashed_filename(x, k, l; prefix="BASISData_", suffix=".h5")

Return an SHA-256 filename derived from the three grid vectors.
Basis directory caches also include parameters, weights, and custom-basis identity.
"""
function generate_hashed_filename(x::AbstractVector, k::AbstractVector, l::AbstractVector;
                                  prefix::String="BASISData_", suffix::String=".h5")
    return hashed_filename((collect(Int, x), collect(k), collect(Int, l)); prefix, suffix)
end

function hashed_filename(values; prefix="BASISData_", suffix=".h5")
    buffer = IOBuffer()
    serialize(buffer, values)
    return prefix * bytes2hex(sha256(take!(buffer))) * suffix
end

end # module PBASEXTypes
