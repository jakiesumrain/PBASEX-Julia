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

function cubic_bspline(r, k, params)
    t = abs((r - k) / params[1])
    if t < 1
        return (4 - 6t^2 + 3t^3) / 6
    elseif t < 2
        return (2 - t)^3 / 6
    end
    return zero(t)
end

ray_distance(r, R) = r > R ? sqrt((r - R) * (r + R)) : zero(R)

cubic_bspline_limit(R, k, params) = ray_distance(k + 2 * params[1], R)

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

function projected_integral(f::GIntegrand{T}, zIP, rtol) where T
    limit = zIP(f.R_pixel, f.k_val, f.params)
    limit > sqrt(eps(T)) || return zero(T)
    return first(quadgk(f, zero(T), T(limit); rtol))
end

function projected_integral(f::GIntegrand{T,typeof(cubic_bspline)},
                            ::typeof(cubic_bspline_limit), rtol) where T
    R, k, h = f.R_pixel, f.k_val, f.params[1]
    R < k + 2h || return zero(T)

    # Knots at or inside R all map to u=0. Keep just the last of those.
    knots = ntuple(i -> ray_distance(k + (i - 3) * h, R), 5)
    start = something(findlast(iszero, knots), 1)
    return first(quadgk(f, knots[start:end]...; rtol))
end

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
    - Integrates with `QuadGK.quadgk`. The built-in cubic B-spline uses its support
      boundaries and polynomial knots as the integration endpoints.
- Stores the integration results in an intermediate matrix.
- Applies optional weights, transposes the matrix, and normalizes by `1 / (2*π)`.
- After computing the intermediate weighted matrix, it calculates the SVD `G = U*diagm(S)*V'` 
  and returns U, S, V directly.
"""
function findG(X::AbstractVector, K::AbstractVector{T}, L::AbstractVector{Int},
             rBF::F, zIP::Z, params::AbstractVector{T},
             weights::Union{AbstractMatrix{T}, Nothing}) where {T<:AbstractFloat,F<:Function,Z<:Function}

    # Get lengths of input vectors
    lenX = length(X)
    lenK = length(K)
    lenL = length(L)
    num_pixels = lenX^2
    num_basis_funcs = lenK * lenL

    # Create coordinate grid for the pixels
    # meshgrid provides X coordinates in each row, Y coordinates in each column
    xl_grid, yl_grid = meshgrid(T.(X), T.(X))
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
    progress = Progress(num_pixels; dt=1, desc="Calculating G basis functions...")
    rtol = max(T(1e-8), eps(T))

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

                integrand = GIntegrand(R_pixel, y_pixel, l_val, k_val, params, rBF, leg)
                G_intermediate[basis_idx, ind] = projected_integral(integrand, zIP, rtol)
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
    factors = svd!(G_intermediate)
    println("SVD complete.")
    
    # Return the SVD components in the order equivalent to svd(transpose(scaled G_intermediate))
    return factors.V, factors.S, factors.U
end

end # module BasisFunctions
