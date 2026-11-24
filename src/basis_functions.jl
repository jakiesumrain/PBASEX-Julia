"""
    BasisFunctions

Project radial and angular basis functions onto the VMI detector and factor
the forward matrix `G`.
"""
module BasisFunctions

using Base.Threads
using LinearAlgebra
using QuadGK
using ProgressMeter

using ..MathUtils: meshgrid, leg

export findG

function cubic_bspline(radius, radial_center, basis_parameters)
    knot_offset = abs((radius - radial_center) / basis_parameters[1])
    if knot_offset < 1
        return (4 - 6knot_offset^2 + 3knot_offset^3) / 6
    elseif knot_offset < 2
        return (2 - knot_offset)^3 / 6
    end
    return zero(knot_offset)
end

# Map a spherical radius to its intersection with the line of sight.
ray_distance(radius, detector_radius) = radius > detector_radius ? sqrt((radius - detector_radius) * (radius + detector_radius)) : zero(detector_radius)

cubic_bspline_limit(detector_radius, radial_center, basis_parameters) =
    ray_distance(radial_center + 2 * basis_parameters[1], detector_radius)

struct GIntegrand{T, RadialProfile, LegendrePolynomial}
    detector_radius::T
    polarization_coordinate::T
    legendre_order::Int
    radial_center::T
    basis_parameters::Vector{T}
    radial_profile::RadialProfile
    legendre_polynomial::LegendrePolynomial
end

function (integrand::GIntegrand)(line_of_sight)
    T = eltype(integrand.basis_parameters)
    radius_squared = line_of_sight^2 + integrand.detector_radius^2
    radius = sqrt(radius_squared)

    # The detector y-axis is the polarization axis. Use cos(θ)=0 at the origin.
    cos_theta = radius_squared > zero(T) ? (integrand.polarization_coordinate / radius) : zero(T)
    angular_value = integrand.legendre_polynomial(integrand.legendre_order, cos_theta)
    radial_value = integrand.radial_profile(radius, integrand.radial_center, integrand.basis_parameters)
    return radial_value * angular_value
end

function projected_integral(integrand::GIntegrand{T}, line_of_sight_limit, rtol) where T
    upper_limit = line_of_sight_limit(integrand.detector_radius, integrand.radial_center, integrand.basis_parameters)
    upper_limit > sqrt(eps(T)) || return zero(T)
    return first(quadgk(integrand, zero(T), T(upper_limit); rtol))
end

function projected_integral(integrand::GIntegrand{T,typeof(cubic_bspline)},
                            ::typeof(cubic_bspline_limit), rtol) where T
    detector_radius = integrand.detector_radius
    radial_center = integrand.radial_center
    knot_spacing = integrand.basis_parameters[1]
    detector_radius < radial_center + 2knot_spacing || return zero(T)

    # A ray begins at u=0. Keep a single endpoint for knots at or inside its detector radius.
    ray_knots = ntuple(i -> ray_distance(radial_center + (i - 3) * knot_spacing, detector_radius), 5)
    first_endpoint = something(findlast(iszero, ray_knots), 1)
    return first(quadgk(integrand, ray_knots[first_endpoint:end]...; rtol))
end

"""
    findG(detector_coordinates, radial_centers, legendre_orders,
          radial_profile, line_of_sight_limit, basis_parameters, weights)

Build the forward operator `G` mapping basis coefficients to a folded detector
projection. Each entry integrates a radial profile times `Pₗ(cos θ)` along
one half of the line of sight and applies the CPBASEX normalization `1/(2π)`.
The built-in cubic B-spline splits the integral at its projected knots.

# Arguments
- `detector_coordinates` gives pixel coordinates along each detector axis.
- `radial_centers` gives the radial basis centers in pixels.
- `legendre_orders` gives the angular basis orders.
- `radial_profile(radius, radial_center, basis_parameters)` evaluates a radial basis function.
- `line_of_sight_limit(detector_radius, radial_center, basis_parameters)` supplies the upper integration limit.
- `basis_parameters` holds radial shape parameters, such as Gaussian width or spline knot spacing.
- `weights` is a matrix on the folded detector grid or `nothing`.

# Returns
The thin SVD `(U, S, V)` of `G`, with `G = U * Diagonal(S) * V'`. Rows of
`G` follow column-major detector order. The radial center varies fastest
among its columns.
"""
function findG(detector_coordinates::AbstractVector, radial_centers::AbstractVector{T},
               legendre_orders::AbstractVector{Int}, radial_profile::F, line_of_sight_limit::Z,
               basis_parameters::AbstractVector{T},
               weights::Union{AbstractMatrix{T}, Nothing}) where {T<:AbstractFloat,F<:Function,Z<:Function}
    detector_side = length(detector_coordinates)
    radial_count = length(radial_centers)
    angular_count = length(legendre_orders)
    pixel_count = detector_side^2
    basis_count = radial_count * angular_count

    detector_x_grid, detector_y_grid = meshgrid(T.(detector_coordinates), T.(detector_coordinates))
    detector_x = vec(detector_x_grid)
    detector_y = vec(detector_y_grid)
    detector_radii = @. sqrt(detector_x^2 + detector_y^2)

    # Store G transposed so each thread writes a contiguous column for its pixel.
    forward_transpose = zeros(T, basis_count, pixel_count)
    progress = Progress(pixel_count; dt=1, desc="Calculating G basis functions...")
    rtol = max(T(1e-8), eps(T))
    Threads.@threads for pixel_index = 1:pixel_count
        detector_radius = detector_radii[pixel_index]
        polarization_coordinate = detector_y[pixel_index]

        for angular_index = 1:angular_count
            legendre_order = legendre_orders[angular_index]
            for radial_index = 1:radial_count
                radial_center = radial_centers[radial_index]
                basis_index = (angular_index - 1) * radial_count + radial_index
                integrand = GIntegrand(detector_radius, polarization_coordinate, legendre_order,
                                       radial_center, basis_parameters, radial_profile, leg)
                forward_transpose[basis_index, pixel_index] = projected_integral(integrand, line_of_sight_limit, rtol)
            end
        end
        next!(progress)
    end
    finish!(progress)

    # In this layout, a detector weight scales a column rather than a row.
    if !isnothing(weights)
        @assert size(weights) == (detector_side, detector_side) "Weights matrix must have size ($detector_side, $detector_side)"
        forward_transpose .*= vec(weights)'
    end

    projection_normalization = (2 * T(π))
    forward_transpose ./= projection_normalization

    println("Performing SVD on scaled forward matrix transpose (in-place)...")
    factors = svd!(forward_transpose)
    println("SVD complete.")

    # Swap left and right singular vectors to return the SVD of G.
    return factors.V, factors.S, factors.U
end

end # module BasisFunctions
