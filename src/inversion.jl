"""
    PBASEX

Provides the main `pbasex` function that implements the core PBASEX (Polar Basis 
Set Expansion) algorithm. This algorithm analyzes 2D projection images (like 
photoelectron or photoion images) to reconstruct the underlying 3D distribution, 
extracting energy spectra and angular anisotropy parameters.

It utilizes helper modules for type definitions (`PBASEXTypes`), image processing 
(`ImageProcessing`), basis function calculations (implicitly via `BASISData`), and 
mathematical utilities.
"""
module PBASEX

using LinearAlgebra

using ..PBASEXTypes: BASISData, PBASEXResult
using ..ImageProcessing: foldQuadrant, resizeFolded, unfoldQuadrant

export pbasex

"""
    pbasex(image::AbstractMatrix{T}; 
           basis::Union{AbstractString, Nothing} = nothing, 
           options...) where T<:AbstractFloat -> PBASEXResult{T}

Perform the pBASEX inversion on a 2D image (e.g., photoelectron image) 
to extract the radial (energy) distribution and angular anisotropy parameters (β).

This function first preprocesses the 2D experimental image (folding, resizing, 
applying weights) using `pre_process_image`. It then uses a precomputed or 
newly calculated basis set (`BASISData`) to solve the inverse problem, finding 
the coefficients that best reconstruct the image. From these coefficients, the 
energy spectrum (`IE`) and anisotropy parameters (`betas`) are calculated.

# Arguments
- `image::AbstractMatrix{T}`: 2D array containing the experimental image. Must be square.
- `basis::Union{AbstractString, Nothing}`: Optional path to load/save `BASISData`.
    - If a file path: Loads basis data from this HDF5 file.
    - If a directory path: Looks for/saves a uniquely named basis file (based on hash of `x, k, l`) within this directory.
    - If `nothing` (default): Uses/creates a `data` subdirectory within the package's `src` directory.
- `options...`: Keyword arguments controlling basis definition and processing. See `pre_process_image` 
  for detailed descriptions of the accepted keywords (`x`, `k`, `l`, `params`, `rBF`, `zIP`, 
  `make_images`, `reg`, `center`, `quadrant_filter`, `weights`, `alpha`).

# Returns
- `PBASEXResult{T}`: A structure containing the results:
    - `E`: Energy scale vector (α * x²).
    - `IE`: Intensity vs. Energy vector (radial distribution, `I(E)`).
    - `betas`: Matrix of anisotropy parameters (`β_l(E)`), where `l` corresponds to orders in `basis_data.l` > 0.
    - `reconstructed_image`: The 2D image reconstructed from the calculated coefficients and basis (`G = U*S*V'`), if `make_images=true`.
    - `inverted_image`: The 2D image reconstructed directly from coefficients and inverse basis (`Ginv`), often representing the central slice (Abel inverted), if `make_images=true`.
    - `residuals`: `image - reconstructed_image`, if `make_images=true`.
    - `coefficients`: Matrix of raw basis coefficients `c_kl`, dimensions `(length(k), length(l))`. 

# Algorithm Steps:
1. Call `pre_process_image` to validate inputs, parse options, handle paths, 
   and perform image/weight preprocessing (folding, resizing, weighting).
2. Load or compute `BASISData` using the path and options returned by `pre_process_image`.
3. Vectorize the preprocessed (weighted, resized) image.
4. Calculate coefficients `ckl` using SVD components from `BASISData` and the regularization parameter.
5. Reshape `ckl` into a matrix `(nk, nl)`.
6. Calculate Energy scale `E = α * x^2`.
7. Calculate Intensity `IE` and Betas using `ckl` and `basis_data.frk`.
8. If `make_images=true`:
    - Reconstruct the projected image using `U, S, V`.
    - Reconstruct the inverted (slice) image using `Ginv`.
    - Calculate residuals by comparing the *original* input `image` with the unfolded reconstructed image.
9. Package results into `PBASEXResult` structure.
"""
function pbasex(image::AbstractMatrix{T}; 
                basis::Union{AbstractString, Nothing} = nothing, 
                options...) where T<:AbstractFloat  
    # --- Pre-Processing --- 
    weighted_resized_image, basis_path, basis_options, processing_options = pre_process_image(image, basis, options)

    # --- Load or Compute BASISData --- 
    # BASISData constructor handles loading from file or computing+saving
    # It uses basis_options if computing is necessary.
    # We pass the type T to ensure basis matrices match image type.
    basis_data = BASISData{T}(basis_path; basis_options...)

    # Extract key parameters from the loaded/computed basis
    x_basis = basis_data.x  # Radial grid points used by the basis
    k_basis = basis_data.k  # Radial basis function centers
    l_basis = basis_data.l  # Angular basis function orders
    nk = length(k_basis)
    nl = length(l_basis)
    nx_basis = length(x_basis) # Size of the basis grid (rmax + 1 essentially)

    # --- Coefficient Calculation --- 
    # Vectorize the processed image
    weighted_vec_image = vec(weighted_resized_image)
    
    reg = processing_options[:reg]
    alpha = processing_options[:alpha]
    make_images = processing_options[:make_images]

    # Calculate expansion coefficients using the regularized pseudo-inverse based on SVD
    # G = U * S * V' => G_reg_inv = V * inv(S^2 + reg*I) * S * U'
    # ckl = G_reg_inv * weighted_vec_image
    # Note: BASISData stores S as a vector, needs Diagonal
    S_vec = basis_data.S
    S_reg_inv = Diagonal(S_vec ./ (S_vec.^2 .+ reg)) # Regularized inverse of singular values

    # Check matrix dimensions conform: U'(nx*nx, ?) * weighted_vec_image(nx*nx)
    # V(nk*nl, ?) * S_reg_inv(?) * intermediate_vector(?)
    # Make sure U is the correct matrix from SVD (often denoted Up in context)
    U_matrix = basis_data.U # Assuming U is (num_pixels, num_singular_values)
    V_matrix = basis_data.V # Assuming V is (num_basis_funcs, num_singular_values)

    # Calculation: ckl = V * [S / (S^2 + reg)] * U' * weighted_image_vector
    ckl_vec = V_matrix * (S_reg_inv * (U_matrix' * weighted_vec_image))

    # Reshape coefficients into a matrix [nk, nl]
    ckl = reshape(ckl_vec, nk, nl)

    # --- Calculate Physical Quantities (Energy, Intensity, Betas) --- 
    # Energy scale: E = alpha * r^2 (using the basis grid x_basis)
    E = alpha .* x_basis.^2

    # Calculate radial distribution (IE) and anisotropy parameters (betas)
    # This involves projecting coefficients back using the radial basis functions (frk)
    # Formula derived from pbasex paper/implementations (needs careful check)
    # frk is likely (nk, nx), ckl is (nk, nl)
    # Need element-wise multiplication with x and scaling
    # Result should be size (nx, nl)
    # Original MATLAB seems to use: CKL_FKR = 1/(2*alpha) .* X .* (frk * CKL);
    frk_matrix = basis_data.frk # Shape (nk, nx)
    ckl_projected = frk_matrix' * ckl # Shape (nx, nl)
    # Scale result
    scaled_projection = (1 / (2 * alpha)) .* x_basis .* ckl_projected # Element-wise multiplication with x_basis

    # Extract IE (l=0 component) and Betas (l>0 components normalized by IE)
    # Assumes the first column (l=l_basis[1]) corresponds to the isotropic component
    IE = scaled_projection[:, 1] 
    # Normalize other columns by the first column to get betas
    # Handle potential division by zero if IE is zero at some points
    betas = scaled_projection[:, 2:end] ./ IE

    # Replace values corresponding to IE == 0 with NaN
    betas[IE .== 0, :] .= T(NaN)

    # --- Optional Image Reconstruction --- 
    reconstructed_image_full = nothing
    inverted_image_full = nothing
    residuals = nothing
    
    if make_images
        # Reconstruct the *projection* image using G = U * S * V'
        # image_recon_vec = G * ckl_vec = U * Diagonal(S) * V' * ckl_vec
        image_recon_vec = U_matrix * (Diagonal(S_vec) * (V_matrix' * ckl_vec))

        # If weights were used, divide the reconstructed image by them to compare with original
        if !isnothing(basis_options[:weights])
            # Avoid division by zero - replace zeros in weights with NaN or a small number? Using NaN for now.
            resized_weights_vec = vec(basis_options[:weights])
            safe_weights = replace(resized_weights_vec, zero(T) => T(NaN))
            image_recon_vec = image_recon_vec ./ safe_weights 
        end
        # Reshape folded quadrant and unfold
        reconstructed_quadrant = reshape(image_recon_vec, nx_basis, nx_basis)
        reconstructed_image_full = unfoldQuadrant(reconstructed_quadrant; center_offset=true) # Assuming default offset is appropriate

        # Reconstruct the *inverted* (slice) image using Ginv
        # slice_recon_vec = Ginv * ckl_vec 
        Ginv_matrix = basis_data.Ginv # Shape (num_pixels, num_basis_funcs)
        slice_recon_vec = Ginv_matrix * ckl_vec
        # Reshape folded quadrant and unfold
        inverted_quadrant = reshape(slice_recon_vec, nx_basis, nx_basis)
        inverted_image_full = unfoldQuadrant(inverted_quadrant; center_offset=true)

        # Calculate residuals (Original Image - Reconstructed Projection Image)
        # Need to ensure reconstructed image is same size as original input image
        # unfoldQuadrant output size depends on nx_basis and center_offset
        # Note: Comparing against the *original* input image 'image'
        if size(reconstructed_image_full) == size(image) 
            residuals = image .- reconstructed_image_full
        else
            @warn "Size mismatch between original image $(size(image)) and reconstructed image $(size(reconstructed_image_full)). Residuals not calculated."
            # Alternatively, crop/pad reconstructed image?
        end
    end

    # --- Result Packaging --- 
    return PBASEXResult{T}(
        E,                    # Energy scale
        IE,                   # Intensity distribution
        betas,                # Anisotropy parameters
        reconstructed_image_full, # Full reconstructed projection image (or nothing)
        inverted_image_full,      # Full reconstructed inverted (slice) image (or nothing)
        residuals,            # Residuals (or nothing)
        ckl                   # Raw coefficients matrix [nk, nl]
    )
end


"""
    pre_process_image(image::AbstractMatrix{T}, basis_file::Union{AbstractString, Nothing}, options::Base.Pairs)
        where T<:AbstractFloat 
        -> Tuple{AbstractMatrix{T}, String, Dict{Symbol, Any}, Dict{Symbol, Any}}

Handles input validation, option parsing, path management, and image/weight 
preprocessing for the main `pbasex` function.

# Arguments
- `image::AbstractMatrix{T}`: 2D input image (must be square).
- `basis_file::Union{AbstractString, Nothing}`: Path for loading/saving basis data (file, directory, or nothing).
- `options::Base.Pairs`: Keyword arguments passed from `pbasex`.
    Known keywords:
        - Basis Definition: `:x`, `:k`, `:l`, `:weights`, `:params`, `:rBF`, `:zIP`
            - Be very careful with the basis definition, it can be tricky! It's your responsibility to ensure consistency between the basis definition and the image size.
            - You cannot reuse the svd of the orignal coefficient matrix to calculate the svd of the weighted coefficient matrix.
        - Processing: `:make_images`, `:reg`, `:center`, `:quadrant_filter`, `:alpha`

# Returns
- `Tuple`: A tuple containing:
    - `weighted_resized_image::AbstractMatrix{T}`: The input image after folding, 
      resizing (to match basis grid `nx`), and applying weights (if provided).
    - `basis_path::String`: The resolved absolute path (directory or file) for BASISData.
    - `basis_options::Dict{Symbol, Any}`: Dictionary containing options relevant 
      for the `BASISData` constructor (x, k, l, resized_weights, params, rBF, zIP).
    - `processing_options::Dict{Symbol, Any}`: Dictionary containing options relevant 
      for the main `pbasex` processing steps (make_images, reg, alpha).

# Processing Steps:
1. Validates image dimensions and keyword arguments.
2. Parses basis definition and processing keywords, setting defaults where needed.
3. Validates processing options (e.g., `reg >= 0`, `center` in bounds).
4. Folds the input `image` based on `quadrant_filter` and `center`.
5. Resizes the folded image to match the basis grid size `nx` (derived from `x` option).
6. If `weights` are provided, folds and resizes them similarly.
7. Applies the resized weights to the resized image element-wise.
8. Resolves the `basis_path` (handling default directory creation).
9. Bundles relevant options into dictionaries for return.
"""
function pre_process_image(image::AbstractMatrix{T}, basis_file::Union{AbstractString, Nothing}, options) where T<:AbstractFloat
    # Validate image size
    rows, cols = size(image)
    @assert rows == cols "Input image must be square, got dimensions ($rows, $cols)"
    img_size = rows

    # --- Define and Parse Options --- 
    # Define known option keys
    basis_keys = Set([:x, :k, :l, :weights, :params, :rBF, :zIP])
    processing_keys = Set([:make_images, :reg, :center, :quadrant_filter, :alpha])
    all_known_keys = union(basis_keys, processing_keys)

    # Check for unknown options
    unknown_keys = setdiff(keys(options), all_known_keys)
    if !isempty(unknown_keys)
        @warn "Unknown keyword arguments provided: $unknown_keys. These will be ignored."
    end

    x = get(options, :x, collect(0:511))::Vector{Int}
    k = get(options, :k, x[2:4:end] .+ 0.5)::Vector{T}
    l = get(options, :l, [0, 2, 4])::Vector{Int}
    # Default params depends on k and x
    default_params_val = T[0.7 * length(x) / length(k)] # Avoid division by zero if k is empty
    params = get(options, :params, default_params_val)::Vector{T}
    # Default functions
    default_rBF = (x_val, k_val, p) -> exp(-(x_val - k_val)^2 / (2 * p[1]^2))
    rBF = get(options, :rBF, default_rBF)::Function
    default_zIP = (x_val, k_val, p) -> sqrt(max(T(0.0), (sqrt(T(20.0))*p[1] + k_val)^2 - x_val^2))
    zIP = get(options, :zIP, default_zIP)::Function
    # Weights are handled later after resizing
    weights_input = get(options, :weights, nothing)::Union{AbstractMatrix{T}, Nothing}

    nx = length(x) # Size of the basis grid (rmax + 1 essentially)

    # --- Extract Processing Options (with defaults) ---
    make_images = get(options, :make_images, true)::Bool
    reg = T(get(options, :reg, zero(T)))::T 
    default_center_x = div(cols, 2) + 1
    default_center_y = div(rows, 2) + 1
    center = get(options, :center, (default_center_x, default_center_y))::Tuple{Int, Int}
    quadrant_filter = get(options, :quadrant_filter, [1, 1, 1, 1])::Vector{Int}
    alpha = T(get(options, :alpha, one(T)))::T
    
    # --- Validate Processing Options --- 
    @assert reg >= zero(T) "Regularization parameter reg=$reg must be non-negative"
    @assert length(quadrant_filter) == 4 && all(val -> val in (0, 1), quadrant_filter) "quadrant_filter must be a 4-element vector of 0s and 1s"
    cx, cy = center
    @assert 1 <= cx <= cols "Center x=$cx is outside image column bounds [1, $cols]"
    @assert 1 <= cy <= rows "Center y=$cy is outside image row bounds [1, $rows]"
    @assert alpha > zero(T) "Energy scaling factor alpha=$alpha must be positive"
    if !isnothing(weights_input)
        @assert size(weights_input) == size(image) "Input weights matrix must match image size $(size(image)), got $(size(weights_input))"
        if any(weights_input .< 0) @warn "Input weights matrix contains negative values." end
    end

    # --- Image Preprocessing --- 
    folded_image = foldQuadrant(image, center[1], center[2]; quadrant_filter = quadrant_filter)
    resized_image = resizeFolded(folded_image, nx) # Resize to basis grid nx

    # --- Weight Preprocessing --- 
    resized_weights = nothing
    weighted_resized_image = resized_image # Start with unwrapped image
    if !isnothing(weights_input)
        folded_weights = foldQuadrant(weights_input, center[1], center[2]; quadrant_filter = quadrant_filter)
        resized_weights = resizeFolded(folded_weights, nx) # Resize weights to basis grid nx
        # Apply weights (element-wise multiplication)
        weighted_resized_image = resized_weights .* resized_image 
    end

    # --- Construct Option Dictionaries for Return --- 
    # Basis options bundle (note: includes *resized* weights)
    basis_options = Dict{Symbol, Any}(
        :x => x,
        :k => k,
        :l => l,
        :weights => resized_weights, # Pass processed weights to BASISData
        :params => params,
        :rBF => rBF,
        :zIP => zIP
    )

    # Processing options bundle
    processing_options = Dict{Symbol, Any}(
        :make_images => make_images,
        :reg => reg,
        :alpha => alpha
        # Center and quadrant_filter are only used here, not needed later
    )

    # --- Basis Data Path Handling (Simplified) --- 
    # Determine the path string to pass to the BASISData constructor.
    # The constructor itself will handle directory vs file and existence checks.
    local basis_path::String # Ensure basis_path is defined in this scope
    if isnothing(basis_file)
        # Default case: Use the 'data' directory in the package root.
        default_basis_dir = joinpath(dirname(@__DIR__), "data")
        if !isdir(default_basis_dir)
            try
                @info "Default basis directory not found. Creating: $default_basis_dir"
                mkpath(default_basis_dir)
            catch e
                error("Failed to create default basis directory: $default_basis_dir. Reason: $(sprint(showerror, e))")
            end
        end
        basis_path = default_basis_dir
        @info "Using default basis data directory: $basis_path"
    else
        # User provided a path (file or directory)
        # Pass it directly to BASISData constructor, which will resolve it.
        basis_path = basis_file
        # --- Removed existence check: !ispath(basis_file) -> error --- 
        # BASISData constructor now handles non-existent file paths correctly.
    end

    return weighted_resized_image, basis_path, basis_options, processing_options
end
    
    
    
end # module PBASEX
