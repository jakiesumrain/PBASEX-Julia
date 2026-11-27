"""
    InverseBasisFunctions

Evaluate the polar basis on the central slice of the reconstructed distribution.
"""
module InverseBasisFunctions

using LinearAlgebra
using ProgressMeter

using ..MathUtils: meshgrid, leg, axial_coordinates

export build_central_slice_basis

"""
    build_central_slice_basis(detector_coordinates, radial_centers, legendre_orders,
                              radial_profile, basis_parameters; geometry=:quadrant)

Build `Ginv`, the matrix that evaluates the fitted distribution in the central
plane at line-of-sight coordinate zero. Each entry is the radial profile times
`Pₗ(y/r)`, with the polarization axis along `y`.

# Arguments
- `detector_coordinates` gives pixel coordinates along each slice axis.
- `radial_centers` gives the radial basis centers in pixels.
- `legendre_orders` gives the angular basis orders.
- `radial_profile(radius, radial_center, basis_parameters)` evaluates a radial basis function.
- `basis_parameters` holds the radial shape parameters.
- `geometry=:half_plane` samples both axial hemispheres and includes the physical normalization `1/(4π)`.

# Returns
A matrix with one row per slice pixel and one column per basis function.
The radial center varies fastest among the columns, as in `build_forward_basis`.
Multiply `Ginv` by the fitted coefficient vector to recover the slice,
then divide by `4π` for the CPBASEX-compatible quadrant matrix.
"""
function build_central_slice_basis(detector_coordinates::AbstractVector, radial_centers::AbstractVector{T},
                                   legendre_orders::AbstractVector{Int}, radial_profile::Function,
                                   basis_parameters::AbstractVector{T}; geometry=:quadrant) where T<:AbstractFloat
    # Flatten in column order so the radial center varies fastest.
    basis_order_grid, basis_center_grid = meshgrid(legendre_orders, radial_centers)
    basis_orders = vec(basis_order_grid)
    basis_centers = vec(basis_center_grid)

    detector_side = length(detector_coordinates)
    signed_axis = axial_coordinates(detector_coordinates, geometry)
    radial_count = length(radial_centers)
    angular_count = length(legendre_orders)
    pixel_count = length(signed_axis) * detector_side
    basis_count = radial_count * angular_count
    @assert length(basis_orders) == basis_count "Meshgrid size mismatch for basis functions"

    slice_x_grid, slice_y_grid = meshgrid(T.(detector_coordinates), T.(signed_axis))
    slice_x = vec(slice_x_grid)
    slice_y = vec(slice_y_grid)
    slice_radii = @. sqrt(slice_x^2 + slice_y^2)
    # At the origin, cos(θ) is undefined. P₀ still returns one.
    cos_theta = slice_y ./ slice_radii

    slice_basis = zeros(T, pixel_count, basis_count)
    progress = Progress(basis_count; dt=1, desc="Calculating G_inv basis functions...")
    Threads.@threads for basis_index = 1:basis_count
        radial_center = basis_centers[basis_index]
        legendre_order = basis_orders[basis_index]
        radial_values = radial_profile.(slice_radii, radial_center, Ref(basis_parameters))
        angular_values = leg(legendre_order, cos_theta)
        @. slice_basis[:, basis_index] = radial_values * angular_values
        next!(progress)
    end
    finish!(progress)

    # Match MATLAB's zero anisotropic slice terms at the undefined origin.
    replace!(slice_basis, T(NaN) => zero(T))
    geometry === :half_plane && (slice_basis ./= 4 * T(π))
    return slice_basis
end

end # module InverseBasisFunctions
