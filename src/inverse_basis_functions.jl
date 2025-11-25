"""
    InverseBasisFunctions

Provides functions for calculating the inverse basis transformation matrix `Ginv`
used in the PBASEX algorithm. This matrix projects the measured 2D image data
onto the basis set to determine the expansion coefficients.
"""
module InverseBasisFunctions

using LinearAlgebra
using ProgressMeter # For displaying calculation progress

# Assuming MathUtils provides meshgrid and leg (Legendre polynomials)
include("math_utils.jl")
using .MathUtils: meshgrid, leg

export findGinv

"""
    findGinv(X::AbstractVector, K::AbstractVector{T}, L::AbstractVector{Int}, 
             rBF::Function, params::AbstractVector{T}) where T<:AbstractFloat -> Matrix{T}

Calculate the inverse basis transformation matrix `Ginv` for the pBASEX algorithm.

This matrix maps from the measured 2D image pixel space back to the pBASEX 
basis coefficients. Each column of `Ginv` represents a basis function 
(a specific combination of radial `k` and angular `l`) evaluated at all 
pixel coordinates.

Multiplying the transpose of this matrix by the flattened image vector yields
the basis coefficients: `coefficients = Ginv' * image_vector` (or equivalent solve).

Unlike the forward matrix `G` (which involves line-of-sight integration), `Ginv`
is computed directly by evaluating the product of the radial basis function (`rBF`)
and the angular basis function (Legendre polynomial `P_l`) at each pixel's 
radius and angle.

# Arguments
- `X::AbstractVector`: Vector of radial pixel coordinates (e.g., `0:rmax`). Defines the image grid.
- `K::AbstractVector{T}`: Vector defining the centers or characteristic radii of the radial basis functions.
- `L::AbstractVector{Int}`: Vector of orders for the Legendre polynomials (angular basis functions).
- `rBF::Function`: The radial basis function, expected signature `rBF(r, k, params)`.
- `params::AbstractVector{T}`: Vector of parameters passed to `rBF`.

# Returns
- `Ginv::Matrix{T}`: The inverse basis transformation matrix. Its dimensions are
  `(num_pixels, num_basis_functions)`, which is `(length(X)^2, length(K) * length(L))`.
  `Ginv[:, j]` contains the j-th basis function evaluated at all pixels.

# Details
- Creates a 2D pixel grid based on `X`.
- Calculates the radius `R` and angle cosine `cosθ = y/R` for each pixel.
- Uses multithreading (`Threads.@threads`) to parallelize the calculation over the basis functions `(k, l)`.
- For each basis function `(k, l)`:
    - Calculates the radial component `rBF.(R, k, params)` for all pixels using broadcasting.
    - Calculates the angular component `leg(l, costh)` for all pixels using broadcasting.
    - Stores the element-wise product `radial .* angular` as a column in the `Ginv` matrix.
- Uses `ProgressMeter` to show progress.
"""
function findGinv(X::AbstractVector, K::AbstractVector{T}, L::AbstractVector{Int}, 
                rBF::Function, params::AbstractVector{T}) where T<:AbstractFloat
    
    # Create a grid of basis function parameters (k, l)
    # lg will have dimensions (lenK, lenL), kg will have dimensions (lenK, lenL)
    # These are flattened later for easier iteration
    lg_grid, kg_grid = meshgrid(L, K)
    lg = vec(lg_grid) # Vector of Legendre orders l
    kg = vec(kg_grid) # Vector of radial basis centers k

    # Get lengths of inputs
    lenX = length(X)
    lenK = length(K)
    lenL = length(L)
    num_pixels = lenX^2
    num_basis_funcs = lenK * lenL
    @assert length(lg) == num_basis_funcs "Meshgrid size mismatch for basis functions"
    
    # --- Calculate Pixel Coordinates --- 
    # Create coordinate grid for the pixels
    xl_grid, yl_grid = meshgrid(X, X)
    # Flatten the grids into vectors
    xl = vec(xl_grid) # Vector of x-coordinates for each pixel
    yl = vec(yl_grid) # Vector of y-coordinates for each pixel
    # Calculate the projected radius R for each pixel
    R = @. sqrt(xl^2 + yl^2)
    # Calculate the cosine of the angle (theta) for each pixel (relative to y-axis)
    costh = yl ./ R
    # Handle R=0 case (origin): cos(theta) is undefined, but for P_l(cosθ), only P_0 is non-zero (and equals 1)
    # We set cos(0) = 0 here; leg(l, 0) will be zero for l!=0, and leg(0, 0) is handled.
    #costh[R .== 0] .= 0.0 
    # Explicitly handle potential NaNs remaining (shouldn't occur with R==0 check)
    #costh[isnan.(costh)] .= 0.0 
    
    # --- Initialize Output Matrix --- 
    # Dimensions: (num_pixels, num_basis_functions)
    # Each column corresponds to one basis function (k, l) evaluated at all pixels.
    Ginv = zeros(T, num_pixels, num_basis_funcs)
    
    # --- Calculate Basis Functions (Parallelized) --- 
    # Setup progress meter for visual feedback
    # Shows progress over the total number of basis functions
    progress = Progress(num_basis_funcs, 1, "Calculating G_inv basis functions...")

    # Use multithreading to compute columns of Ginv in parallel.
    # Each thread handles a subset of the basis functions.
    Threads.@threads for basis_idx = 1:num_basis_funcs
        k = kg[basis_idx] # Radial basis center for this column
        l = lg[basis_idx] # Legendre order for this column

        # Calculate radial component for ALL pixels for this 'k'
        # Result is a vector of size num_pixels.
        radial_component = rBF.(R, k, Ref(params))

        # Calculate angular component for ALL pixels for this 'l'
        # Result is a vector of size num_pixels.
        angular_component = leg(l, costh) # Uses MathUtils.leg

        # Combine components: Ginv column is the product of radial and angular parts
        @. Ginv[:, basis_idx] = radial_component * angular_component

        # Update progress meter (thread-safe)
        next!(progress)
    end

    # Finish progress meter
    finish!(progress)
    
    # Replace any potential NaN values with zeros (as a final safety check)
    Ginv[isnan.(Ginv)] .= 0.0
    
    # Note: The matrix is NOT transposed here.
    # It remains (num_pixels, num_basis_functions).
        
    return Ginv
end

end # module InverseBasisFunctions