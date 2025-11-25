module PBASEXTypes

using HDF5
using Serialization
using SHA
using LinearAlgebra # Needed for svd

# Include necessary functions (assuming they are defined elsewhere, e.g., basis_functions.jl)
# If these are defined within this module or another included module, adjust accordingly.
# We assume findG and findGinv are available in the scope where BASISData constructors are called.
# If they are in other modules, they might need qualification (e.g., BasisFunctions.findG)
# For now, we'll assume they are accessible.
include("basis_functions.jl") # Assuming findG is here
include("inverse_basis_functions.jl") # Assuming findGinv is here
include("math_utils.jl") # Assuming meshgrid is here
using .BasisFunctions: findG # Explicitly import if needed
using .InverseBasisFunctions: findGinv # Explicitly import if needed
using .MathUtils: meshgrid # Explicitly import if needed

export PBASEXResult, BASISData, generate_hashed_filename

"""
    PBASEXResult{T<:AbstractFloat}

Stores the results obtained from applying the pBASEX algorithm to one or more images.

This structure encapsulates the primary outputs, including the energy distribution,
anisotropy parameters, expansion coefficients, and reconstructed/inverted images.

# Fields
- `E::Vector{T}`: Energy scale corresponding to the rows of `IE` and `betas`. Often proportional to radius squared (r²).
- `IE::Matrix{T}`: Intensity Energy spectrum. `IE[i, j]` is the intensity at energy `E[i]` for the j-th input image.
- `betas::Matrix{T}`: Angular anisotropy parameters (β values). `betas[i, k]` contains the β_k parameter at energy `E[i]`. Note: The index `k` depends on the Legendre polynomials `l` used in the basis.
- `reconstructed_image::Union{Matrix{T}, Nothing}`: The image reconstructed from the calculated coefficients and basis functions. `nothing` if reconstruction was skipped.
- `inverted_image::Union{Matrix{T}, Nothing}`: The central-slice image obtained from the Abel inversion of the reconstructed image (if applicable). `nothing` if inversion was skipped or not applicable.
- `residuals::Union{Matrix{T}, Nothing}`: Difference between the original processed image and the `reconstructed_image`. `nothing` if reconstruction was skipped.
- `coefficients::Matrix{T}`: Raw expansion coefficients obtained from the fit. `coefficients[i, j]` corresponds to the i-th basis function coefficient for the j-th input image.

# Notes
- The dimensions of `IE`, `betas`, and `coefficients` depend on the number of energy points (derived from the basis `k` vector) and the number of input images processed.
- The interpretation of `betas` depends on the angular basis functions (Legendre polynomials `l`) used. For a typical basis `l = [0, 2, 4]`, `betas[:, 1]` would correspond to β₂, and `betas[:, 2]` to β₄.
"""
struct PBASEXResult{T<:AbstractFloat}
    E::Vector{T}
    IE::Vector{T}
    betas::Matrix{T} # Note: Dimension corrected from Array{T,2} to Matrix{T}
    reconstructed_image::Union{Matrix{T}, Nothing}
    inverted_image::Union{Matrix{T}, Nothing}
    residuals::Union{Matrix{T}, Nothing}
    coefficients::Matrix{T}
end

"""
    BASISData{T<:AbstractFloat}

Represents the pBASEX basis set and associated precomputed matrices.

This structure holds the definition of the radial (`k`) and angular (`l`) basis functions
over the image coordinate grid (`x`). It also stores potentially large matrices derived
from these definitions (like the forward transform `G`, its SVD components `U, S, V`,
the inverse transform `Ginv`, and the radial basis functions evaluated on the grid `frk`)
to avoid redundant calculations.

Instances can be created by:
1. Providing all parameters and letting the constructor compute the matrices.
2. Providing parameters and a file path to save the computed matrices.
3. Providing a file path to load precomputed matrices.

# Fields
- `x::Vector{Int}`: Grid points for the image radius (or relevant coordinate).
- `k::Vector{T}`: Centers (or characteristic points) of the radial basis functions.
- `l::Vector{Int}`: Orders of the Legendre polynomials used for the angular basis.
- `weights::Union{Matrix{T}, Nothing}`: Optional weighting matrix applied during basis generation.
- `params::Vector{T}`: Parameters controlling the shape/width of the radial basis functions (used by `rBF`).
- `rBF::Function`: The radial basis function, `rBF(x, k, params)`.
- `zIP::Function`: Function related to the Abel transform kernel, `zIP(x, k, params)`.
- `U::Matrix{T}`: Left singular vectors of the forward transform matrix `G`.
- `S::Vector{T}`: Singular values of `G`. Stored as a vector.
- `V::Matrix{T}`: Right singular vectors of `G` (as `V`, not `V'`).
- `frk::Matrix{T}`: Radial basis functions evaluated on the `x, k` grid: `frk[i, j] = rBF(x[i], k[j], params)`.
- `Ginv::Matrix{T}`: The inverse transformation matrix.
"""
struct BASISData{T<:AbstractFloat}
    # Basis Definition
    x::Vector{Int}
    k::Vector{T}
    l::Vector{Int}

    # Optional Weights
    weights::Union{Matrix{T}, Nothing}

    # Radial Basis Function Parameters and Functions
    params::Vector{T}
    rBF::Function # Radial Basis Function signature: rBF(x, k, params)
    zIP::Function # the integration upper limit

    # Precomputed Matrices (from G = U * Diagonal(S) * V')
    U::Matrix{T}
    S::Vector{T} # Store singular values as a vector
    V::Matrix{T}
    frk::Matrix{T} # Cache of rBF(x, k, params)
    Ginv::Matrix{T} # Inverse transform matrix
end

# --- Constructors for BASISData ---

"""
    BASISData{T}(; kwargs...) where T<:AbstractFloat

Internal parameterized constructor for `BASISData`.

Calculates the necessary matrices (`U`, `S`, `V`, `Ginv`, `frk`) based on the provided
or default basis parameters (`x`, `k`, `l`, `params`, `rBF`, `zIP`, `weights`).
This constructor is typically called by the user-facing constructors and focuses 
*only* on computation, not saving.

# Keyword Arguments
- `x`: Vector{Int}, defaults to `0:511`.
- `k`: Vector{T}, defaults based on `x` (`x[2:4:end] .+ 0.5`).
- `l`: Vector{Int}, defaults to `[0, 2, 4]`.
- `params`: Vector{T}, defaults based on `k` and `x` (`[0.7 * length(k) / length(x)]`).
- `weights`: Matrix{T} or Nothing, defaults to `nothing`.
- `rBF`: Function, defaults to a Gaussian: `(x, k, p) -> exp(-(x - k)^2 / (2 * p[1]^2))`.
- `zIP`: Function, defaults to `(x, k, p) -> sqrt((sqrt(2*10)*p[1] + k)^2 - x^2)`.
- `U`, `S`, `V`: Optional precomputed SVD components of `G`. If not provided, `findG` is called.
- `Ginv`: Optional precomputed inverse transform matrix. If not provided, `findGinv` is called.
- `frk`: Optional precomputed radial basis function cache. If not provided, it's computed.
# Note: `basis_file` keyword is NO LONGER used here. Saving is handled by the other constructor.
"""
function BASISData{T}(;kwargs...) where T<:AbstractFloat

    # --- Process Keyword Arguments and Set Defaults ---
    x = get(kwargs, :x, collect(0:511))::Vector{Int}
    l = get(kwargs, :l, [0, 2, 4])::Vector{Int}

    # Default 'k' depends on 'x'
    k = get(kwargs, :k, T.(x[2:4:end] .+ 0.5))::Vector{T}
    # Default 'params' depends on 'k' and 'x'
    params = get(kwargs, :params, T[0.7 * length(x) / length(k)])::Vector{T}

    # Optional fields
    weights = get(kwargs, :weights, nothing)::Union{Matrix{T}, Nothing}
    # Default Radial Basis Function (Gaussian)
    rBF = get(kwargs, :rBF, (x, k, prms) -> exp(-(x - k)^2 / (2 * prms[1]^2)))::Function
    # Default Inverse Problem related function
    zIP = get(kwargs, :zIP, (x, k, prms) -> sqrt(max(T(0.0), (sqrt(T(20.0))*prms[1] + k)^2 - x^2)))::Function # Ensure type stability

    local U::Matrix{T}, S::Vector{T}, V::Matrix{T}, frk::Matrix{T}, Ginv::Matrix{T}

    # --- Calculate or Retrieve Matrices ---
    # SVD is now done within findG
    if !haskey(kwargs, :U) || !haskey(kwargs, :S) || !haskey(kwargs, :V)
        # Call the modified findG which returns U, S, V directly
        println("Calling findG to calculate basis and perform SVD...") # Add feedback
        U, S, V = findG(x, k, l, rBF, zIP, params, weights)
        println("findG calculation complete.")
        # # --- SVD calculation removed from here --- 
        # G = findG(x, k, l, rBF, zIP, params, weights) 
        # F = svd(G) 
        # U = F.U 
        # S = F.S 
        # V = F.V 
    else
        # Keys :U, :S, :V are guaranteed to exist here, use direct indexing.
        U = kwargs[:U]::Matrix{T}
        S = kwargs[:S]::Vector{T} # Expect Vector
        V = kwargs[:V]::Matrix{T}
    end

    # Calculate Ginv if not provided
    if !haskey(kwargs, :Ginv)
        # Calculate the inverse transform matrix
        # Assumes findGinv is accessible here
        Ginv = findGinv(x, k, l, rBF, params)
    else
        # Key :Ginv is guaranteed to exist here, use direct indexing.
        Ginv = kwargs[:Ginv]::Matrix{T}
    end

    # Calculate frk (cached rBF evaluations) if not provided
    if !haskey(kwargs, :frk)
        # Assumes meshgrid is accessible here
        R, K = meshgrid(T.(x), k) # Ensure type consistency
        frk = rBF.(R, K, Ref(params)) # Use broadcasting with Ref for params
    else
        # Key :frk is guaranteed to exist here, use direct indexing.
        frk = kwargs[:frk]::Matrix{T}
    end

    # --- Create and Return Instance ---
    return BASISData{T}(
        x, k, l, weights, params,
        rBF, zIP, U, S, V, frk, Ginv
    )
end

"""
    BASISData{T}(basis_path::AbstractString; kwargs...) where T<:AbstractFloat
    BASISData(basis_path::AbstractString; kwargs...)

User-facing constructor for `BASISData` that handles loading from or saving to a file.

Determines whether to load existing basis data or compute new data based on the
provided `basis_path` and basis parameters.

# Arguments
- `basis_path::AbstractString`: Path to the basis data file or directory.
  - If `basis_path` is an existing HDF5 file: Attempts to load the `BASISData` from it.
    Ignores other keyword arguments defining basis parameters (`x`, `k`, `l`, etc.)
    as the loaded file defines the basis.
  - If `basis_path` is an existing directory: Constructs a unique filename within that
    directory based on the hash of `x`, `k`, `l` (provided via `kwargs` or defaults).
    - If the file exists, loads it.
    - If the file does not exist, computes the basis data using `kwargs` (or defaults)
      by calling the computation constructor `BASISData{T}(;kwargs...)`, and then 
      attempts to save it to the generated filename using `save_basis_data`.
  - If `basis_path` is neither an existing file nor directory, throws an `ArgumentError`.

- `kwargs...`: Keyword arguments defining the basis (`x`, `k`, `l`, `params`, `rBF`, `zIP`, `weights`).
  These are used **only** when computing a new basis (i.e., when `basis_path` is a directory
  and the hashed file doesn't exist).

# Returns
- A `BASISData{T}` object, either loaded from file or newly computed (and potentially saved).

# Notes
- The non-parameterized version `BASISData(basis_path; kwargs...)` defaults `T` to `Float64`.
- Relies on `generate_hashed_filename` to create consistent filenames based on `x`, `k`, `l`.
- Loading assumes the HDF5 file contains the required fields (`x`, `k`, `l`, `U`, `S`, `V`, `frk`, `Ginv`).
- Saving uses the `save_basis_data` function (assumed to be defined elsewhere in this module).
"""
function BASISData{T}(basis_path::AbstractString; kwargs...) where T<:AbstractFloat

    # Define allowed keywords for creating a *new* basis, used for validation
    creation_keys = Set([:x, :k, :l, :weights, :params, :rBF, :zIP])
    unknown_keys = setdiff(keys(kwargs), creation_keys)
    if !isempty(unknown_keys)
        @warn "Unknown keyword arguments provided to BASISData file constructor (will be ignored if loading): $unknown_keys"
    end

    fullfilename = ""
    load_attempt = false # Flag to indicate if we should try loading

    if isdir(basis_path)
        # Case 1: User provided a directory. Generate hashed filename for caching.
        @info "Basis path provided is a directory. Using hash-based caching..."
        # Get x, k, l from kwargs or defaults to generate the hash
        temp_x = get(kwargs, :x, collect(0:511))
        temp_l = get(kwargs, :l, [0, 2, 4])
        # Ensure k is of type T for hash consistency if not provided
        temp_k_maybe_T = get(kwargs, :k, T.(temp_x[2:4:end] .+ 0.5))
        # Use generate_hashed_filename (ensure it uses T for k)
        filename = generate_hashed_filename(temp_x, temp_k_maybe_T, temp_l)
        fullfilename = joinpath(basis_path, filename)
        load_attempt = isfile(fullfilename) # Try loading only if the specific hashed file exists
        # Proceed to loading/computation block
    else
        # Case 2: User provided path is not a directory. Assume it's a specific file path.
        fullfilename = basis_path 
        # Check if parent directory exists, otherwise saving will fail
        parent_dir = dirname(fullfilename)
        if !isdir(parent_dir)
            # Check if parent_dir is empty string (meaning relative path in current dir)
            # If empty, the current directory must exist, which is implicitly true. 
            # Otherwise, throw error.
            if parent_dir != ""
                throw(ArgumentError("Parent directory for specified basis_file does not exist: $(parent_dir)"))
            end
        end
        
        load_attempt = isfile(fullfilename)
        
        if load_attempt && !isempty(kwargs)
             # Warn only if loading an existing specified file and kwargs were also given
            @info "Basis parameters provided via keyword arguments will be ignored when loading directly from specified file: $(fullfilename)"
        elseif !load_attempt
            @info "Specified basis_file does not exist: $(fullfilename). Will compute and save to this path."
        end
        # Proceed to loading/computation block
    end
    
    # --- Attempt Loading or Compute/Save ---
    if load_attempt
        # Try loading from the determined fullfilename
        try
            @info "Attempting to load BASISData from file: $fullfilename"
            # Initialize variables to be loaded
            local x_loaded, k_loaded, l_loaded, U_loaded, S_loaded, V_loaded, frk_loaded, Ginv_loaded
            local weights_loaded = nothing # Assume weights might not be saved
            local params_loaded, rBF_loaded, zIP_loaded # Cannot reliably load functions/params yet

            h5open(fullfilename, "r") do file
                x_loaded = load_required_field(file, "x")
                k_loaded = load_required_field(file, "k")
                l_loaded = load_required_field(file, "l")
                U_loaded = load_required_field(file, "U")
                S_loaded = load_required_field(file, "S") # Should be Vector
                V_loaded = load_required_field(file, "V")
                frk_loaded = load_required_field(file, "frk")
                Ginv_loaded = load_required_field(file, "Ginv")
                # Optionally load weights if they were saved
                # if haskey(file, "weights_present") && read(file["weights_present"])
                #     weights_loaded = load_required_field(file, "weights")
                # end
            end

            # Consistency Checks (Optional but Recommended)
            # Check if loaded types match expected T
            # Check if loaded x, k, l match kwargs if directory mode was used? (Tricky)
            @assert eltype(k_loaded) == T "Loaded 'k' vector type ($(eltype(k_loaded))) does not match expected type $T"
            # Add more assertions for other matrices/vectors if needed

            # Need to provide default/placeholder functions for rBF, zIP and params, weights
            # These are not used if only loading, but the struct requires them.
            # Get defaults matching the loaded x, k, l structure
            default_params = T[0.7 * length(x_loaded) / length(k_loaded)]
            default_rBF = (x, k, p) -> exp(-(x - k)^2 / (2 * p[1]^2))
            default_zIP = (x, k, p) -> sqrt(max(T(0.0), (sqrt(T(20.0))*p[1] + k)^2 - x^2))

            @info "Successfully loaded BASISData from file: $fullfilename"
            # Create BASISData instance with loaded data
            return BASISData{T}(
                x_loaded, k_loaded, l_loaded, weights_loaded, default_params,
                default_rBF, default_zIP, # Use defaults/placeholders
                U_loaded, S_loaded, V_loaded, frk_loaded, Ginv_loaded
            )
        catch e
            @error "Failed to load BASISData from file: $fullfilename. Check file integrity and required fields." exception=(e, catch_backtrace())
            rethrow(e) # Rethrow after logging
        end
    else
        # Case: File does not exist (and path was a directory), compute and save.
        @info "Basis file not found: $fullfilename. Computing new BASISData..."
        # 1. Call the main computation constructor (which no longer saves)
        # It will use kwargs or defaults for computation.
        computed_basis = BASISData{T}(; kwargs...)
        
        # 2. Call the newly added save function
        try
            save_basis_data(computed_basis, fullfilename)
            @info "Computed BASISData saved to file: $fullfilename"
        catch e
            @error "Failed to save newly computed BASISData to file: $fullfilename" exception=(e, catch_backtrace())
            # Decide whether to rethrow, warn, or just proceed returning the computed data
            # rethrow(e) # Option: Stop execution if saving fails
        end
        
        # 3. Return the computed basis object
        return computed_basis 
    end
end

# Non-parameterized version calling the parameterized one with T=Float64
function BASISData(basis_path::AbstractString; kwargs...)
    return BASISData{Float64}(basis_path; kwargs...)
end

"""
    save_basis_data(basis_data::BASISData{T}, basis_file_path::String) where T

Saves the contents of a `BASISData` object to a specified HDF5 file.

This function is typically called by the file-handling `BASISData` constructor 
when a new basis set needs to be computed and cached.

# Arguments
- `basis_data::BASISData{T}`: The computed basis data object to save.
- `basis_file_path::String`: The full path to the HDF5 file where the data should be saved.

# Details
- Opens the `basis_file_path` in write mode (`"w"`), overwriting if it exists.
- Writes the core basis definition fields (`x`, `k`, `l`).
- Writes the computed matrices (`U`, `S`, `V`, `frk`, `Ginv`).
- Optionally could be extended to save `weights` and `params` if needed for 
  consistency checks during loading (currently commented out in typical implementations).
- Uses `h5open` for file handling.

# Errors
- Rethrows any exceptions encountered during file opening or writing (e.g., 
  permissions issues, invalid path, HDF5 library errors).
"""
function save_basis_data(basis_data::BASISData{T}, basis_file_path::String) where T
    try
        h5open(basis_file_path, "w") do file
            # Write core definition
            write(file, "x", basis_data.x)
            write(file, "k", basis_data.k)
            write(file, "l", basis_data.l)
            
            # --- Save weights if they exist --- 
            #if !isnothing(basis_data.weights)
            #     write(file, "weights_present", true)
            #     write(file, "weights", basis_data.weights)
            #else
            #     write(file, "weights_present", false)
            #end
            # Optionally save params too
            # write(file, "params", basis_data.params)
            # --- End Optional --- 

            # Write computed/loaded matrices
            write(file, "U", basis_data.U)
            write(file, "S", basis_data.S) # Vector
            write(file, "V", basis_data.V)
            write(file, "frk", basis_data.frk)
            write(file, "Ginv", basis_data.Ginv)
        end
    catch e
        # Let the caller handle logging/rethrowing based on context
        rethrow(e) 
    end
end

"""
    load_required_field(file::HDF5.File, field::String)

Helper function to read a dataset `field` from an opened HDF5 file `file`.
Throws an `ArgumentError` if the field does not exist.

# Arguments
- `file`: An opened HDF5 file object (e.g., from `h5open`).
- `field`: The name of the dataset to read.

# Returns
- The data read from the specified dataset.
"""
function load_required_field(file::Union{HDF5.File, HDF5.Group}, field::String)
    if !haskey(file, field)
        throw(ArgumentError("Required field '$field' missing from HDF5 file/group: $(HDF5.name(file))"))
    end
    data = read(file[field])
    return data
end


# --- Function to generate the unique filename ---
"""
    generate_hashed_filename(x::AbstractVector, k::AbstractVector, l::AbstractVector;
                              prefix::String="BASISData_", suffix::String=".h5") -> String

Generates a unique, filename-safe string based on the content of vectors x, k, and l
using SHA-256 hashing of their serialized representation.

This ensures that the same basis definition (`x`, `k`, `l`) always results in the same
filename, allowing for reliable caching and retrieval of precomputed basis data.

# Arguments
- `x`, `k`, `l`: The input vectors defining the basis dimensions.
- `prefix`: Optional string to prepend to the hash (default: "BASISData_").
- `suffix`: Optional string to append to the hash (default: ".h5").

# Returns
- A string suitable for use as a filename, uniquely representing the input vectors.
"""
function generate_hashed_filename(x::AbstractVector, k::AbstractVector{T}, l::AbstractVector;
                                  prefix::String="BASISData_", suffix::String=".h5") where T
    # 1. Combine vectors into a single object (Tuple preserves order and types)
    # Ensure types are concrete for consistent serialization
    data_tuple = (Vector{Int}(x), Vector{T}(k), Vector{Int}(l))

    # 2. Serialize the data to an in-memory buffer
    io = IOBuffer()
    Serialization.serialize(io, data_tuple)
    seekstart(io) # Reset buffer position
    serialized_data = take!(io) # Get the bytes
    close(io)

    # 3. Compute the SHA-256 hash of the serialized bytes
    hash_bytes = sha256(serialized_data)

    # 4. Convert the hash bytes to a hexadecimal string
    hash_hex = bytes2hex(hash_bytes)

    # 5. Construct the final filename string
    return prefix * hash_hex * suffix
end

end # module PBASEXTypes