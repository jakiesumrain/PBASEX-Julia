"""
    BasisFunctions

Provides functions for calculating the basis transformation matrices used in the
PBASEX (Polar Basis Set Expansion) algorithm. Specifically, it includes the `findG`
function to compute the forward transformation matrix `G`.
"""
module BasisFunctions

using Base.Threads
using LinearAlgebra
using QuadGK # For numerical integration
using ProgressMeter # For displaying calculation progress

using ..MathUtils: meshgrid, leg

export findG

# --- Define a callable struct for the integrand ---
struct GIntegrand{T, F_rBF, F_leg}
    R_pixel::T
    y_pixel::T
    l_val::Int
    k_val::T
    params::Vector{T}
    rBF::F_rBF   # Store the actual function
    leg::F_leg   # Store the actual function
end

# Define the call method for the struct
function (integrand::GIntegrand)(u)
    T = eltype(integrand.params) # Get the type
    r_3d_squared = u^2 + integrand.R_pixel^2
    # Avoid sqrt(0) which can be slow or problematic
    #if r_3d_squared <= zero(T)
         # Handle degenerate case: contribution is likely zero unless rBF behaves specially at r=0
         # If rBF(0, ...) * leg(...) is non-zero, this needs more care. Assume zero for now.
         #return zero(T)
    #end
    r_3d = sqrt(r_3d_squared)

    # Stable computation for cosine term
    cos_term = r_3d_squared > zero(T) ? (integrand.y_pixel / r_3d) : zero(T)

    angular_part = integrand.leg(integrand.l_val, cos_term)
    radial_part = integrand.rBF(r_3d, integrand.k_val, integrand.params)

    return radial_part * angular_part
end
# --- End of Integrand Struct ---

"""
    findG(X::AbstractVector, K::AbstractVector{T}, L::AbstractVector{Int}, 
          rBF::Function, zIP::Function, params::AbstractVector{T}, 
          weights::Union{AbstractMatrix{T}, Nothing}) where T<:AbstractFloat -> Tuple{Matrix{T}, Vector{T}, Matrix{T}}

Calculate the forward transformation matrix `G` for the pBASEX algorithm 
and return its Singular Value Decomposition (SVD) components (U, S, V).

This matrix maps the pBASEX coefficients (representing the 3D distribution decomposed
into radial and angular basis functions) to the pixel values of the projected 2D image.
It computes the contribution of each basis function (combination of radial `rBF` centered
at `k` and angular `LegendreP(l, cosθ)`) to each pixel `(x, y)` by integrating along
the line-of-sight (related to the Abel transform).

# Arguments
- `X::AbstractVector`: Vector of radial pixel coordinates (e.g., `0:rmax`).
- `K::AbstractVector{T}`: Vector defining the centers or characteristic radii of the radial basis functions.
- `L::AbstractVector{Int}`: Vector of orders for the Legendre polynomials (angular basis functions).
- `rBF::Function`: The radial basis function, expected signature `rBF(r, k, params)` where `r` is radius, `k` is a basis center from `K`, and `params` are shape parameters.
- `zIP::Function`: Function defining the upper limit for the line-of-sight integration, expected signature `zIP(R_pixel, k, params)` where `R_pixel` is the projected radius of the pixel.
- `params::AbstractVector{T}`: Vector of parameters passed to `rBF` and `zIP`.
- `weights::Union{AbstractMatrix{T}, Nothing}`: Optional weighting matrix of the same size as the original image (`len(X) x len(X)`). If provided, its flattened version is applied to the columns of the intermediate `G` matrix.

# Returns
- `Tuple{Matrix{T}, Vector{T}, Matrix{T}}`: A tuple containing:
  - `U::Matrix{T}`: Left singular vectors of `G`.
  - `S::Vector{T}`: Singular values of `G` (as a vector).
  - `V::Matrix{T}`: Right singular vectors of `G` (as `V`, not `V'`).

# Details
- Creates a 2D pixel grid based on `X`.
- For each pixel `(x, y)` and each basis function `(k, l)`:
    - Defines an integrand: `rBF(sqrt(u^2 + R^2), k, params) * P_l(y / sqrt(u^2 + R^2))`
      where `R = sqrt(x^2 + y^2)` is the projected radius and `u` is the integration variable (along the line of sight).
    - Numerically integrates the integrand from `u = 0` to `u = zIP(R, k, params)` using `QuadGK.quadgk`.
- Stores the integration results in an intermediate matrix.
- Applies optional weights, transposes the matrix, and normalizes by `1 / (2*π)`.
- After computing the intermediate weighted matrix, it calculates the SVD `G = U*diagm(S)*V'` 
  and returns U, S, V directly.
"""
function findG(X::AbstractVector, K::AbstractVector{T}, L::AbstractVector{Int}, 
             rBF::Function, zIP::Function, params::AbstractVector{T}, 
             weights::Union{AbstractMatrix{T}, Nothing}) where T<:AbstractFloat

    # Get lengths of input vectors
    lenX = length(X)
    lenK = length(K)
    lenL = length(L)
    num_pixels = lenX^2
    num_basis_funcs = lenK * lenL

    # Create coordinate grid for the pixels
    # meshgrid provides X coordinates in each row, Y coordinates in each column
    xl_grid, yl_grid = meshgrid(X, X) 
    # Flatten the grids into vectors for easier iteration
    xl = vec(xl_grid) # Vector of x-coordinates for each pixel
    yl = vec(yl_grid) # Vector of y-coordinates for each pixel
    # Calculate the projected radius R for each pixel
    R = @. sqrt(xl^2 + yl^2)

    # Initialize the intermediate G matrix.
    # Dimensions: (num_basis_functions, num_pixels)
    # G_intermediate[basis_idx, pixel_idx] will store the integral value.
    G_intermediate = zeros(T, num_basis_funcs, num_pixels)

    # Setup progress meter for visual feedback
    # Shows progress over the total number of basis functions
    progress = Progress(num_pixels, 1, "Calculating G basis functions...")

    # --- Calculate integrals for each pixel and basis function --- 
    Threads.@threads for ind = 1:num_pixels
        R_pixel = R[ind]
        y_pixel = yl[ind]
        
        for lind = 1:lenL
            l_val = L[lind]
            for q = 1:lenK
                k_val = K[q]
                # --- Calculate basis index directly ---
                basis_idx = (lind - 1) * lenK + q

                z_limit_val = zIP(R_pixel, k_val, params)

                # Check if integration is needed (e.g., z_limit might be negative/NaN from custom zIP)
                # Also handle cases where z_limit is very small or zero.
                if !isnan(z_limit_val) && z_limit_val > sqrt(eps(T)) # Integrate only if limit is valid and non-trivial
                    
                    # --- Create an instance of the integrand functor ---
                    integrand_functor = GIntegrand(R_pixel, T(y_pixel), l_val, k_val, params, rBF, leg)

                    # Perform numerical integration
                    result, err = quadgk(integrand_functor, zero(T), z_limit_val, rtol=1e-8)
                    #result, err = quadgk(integrand_functor, zero(T), T(Inf), rtol=1e-8)

                    # Optional: Add back the warning for high integration error if needed
                    # if err > abs(result * 1e-4) ... end

                    G_intermediate[basis_idx, ind] = result
                else
                    # If z_limit is invalid or effectively zero, the integral is zero
                    G_intermediate[basis_idx, ind] = zero(T)
                end
            end # end loop over k (radial basis)
        end # end loop over l (angular basis)
        # Update progress meter (thread-safe)
        next!(progress)
    end # end loop over pixels

    # Finish progress meter
    finish!(progress)
    
    # --- Final processing --- 
    
    # Apply pixel weights if provided (weights are flattened and applied to columns)
    if !isnothing(weights)
        @assert size(weights) == (lenX, lenX) "Weights matrix must have size ($lenX, $lenX)"
        G_intermediate .*= vec(weights)' # Multiply each column by corresponding weight
    end

    # --- Scale G_intermediate in-place and perform SVD --- 
    c = (2 * T(π))
    G_intermediate ./= c # Apply scaling factor in-place
    
    # Perform Singular Value Decomposition on the now-scaled G_intermediate
    # Use svd! for potential memory optimization, as G_intermediate is overwritten
    println("Performing SVD on scaled G_intermediate matrix (in-place)...") # Add feedback
    F = svd!(G_intermediate)
    println("SVD complete.")
    
    # Return the SVD components in the order equivalent to svd(transpose(scaled G_intermediate))
    # U_of_G = F.V
    # S_of_G = F.S
    # V_of_G = F.U
    return F.V, F.S, F.U 
end

end # module BasisFunctions
