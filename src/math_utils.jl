module MathUtils

export leg, meshgrid, cart2pol, pol2cart

"""
    leg(l::Int, x::AbstractArray)

Evaluate Legendre polynomial of order `l` element-wise for an array `x`.

Uses direct formulas for orders l=0 to 4 for efficiency and numerical stability,
and Bonnet's recurrence relation for higher orders.

# Arguments
- `l::Int`: Order of the Legendre polynomial (must be non-negative).
- `x::AbstractArray`: Array of points at which to evaluate the polynomial. Values are typically expected in the range [-1, 1].

# Returns
- An array of the same size as `x` containing the polynomial values P_l(x).

# Examples
```julia
x = [-1.0, 0.0, 1.0]
leg(2, x)  # Returns [1.0, -0.5, 1.0]
```
"""
function leg(order::Int, x::AbstractArray)
    @assert order >= 0 "Order l must be non-negative"
    T = eltype(x)

    if order == 0
        return ones(T, size(x))
    elseif order == 1
        return T.(copy(x))
    elseif order == 2
        return @. (3 * x^2 - 1) / 2
    elseif order == 3
        return @. (5 * x^3 - 3 * x) / 2
    elseif order == 4
        return @. (35 * x^4 - 30 * x^2 + 3) / 8
    end

    previous_legendre = ones(T, size(x))
    current_legendre = T.(copy(x))
    next_legendre = similar(x, T)

    # Carry Pₙ and Pₙ₋₁ through Bonnet’s recurrence.
    for degree in 1:(order-1)
        next_legendre .= ((2*degree+1) .* x .* current_legendre .- degree .* previous_legendre) ./ (degree+1)
        previous_legendre, current_legendre, next_legendre = current_legendre, next_legendre, previous_legendre
    end

    return current_legendre
end

"""
    leg(l::Int, x::Number)

Evaluate Legendre polynomial of order `l` at a single point `x`.

Uses direct formulas for orders l=0 to 4 and Bonnet's recurrence relation
for higher orders.

# Arguments
- `l::Int`: Order of the Legendre polynomial (must be non-negative).
- `x::Number`: Point at which to evaluate the polynomial (typically in range [-1, 1]).

# Returns
- The value of the Legendre polynomial P_l(x).

# Implementation
Uses Bonnet's recurrence relation for l > 4:
(n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)

# Examples
```julia
leg(0, 0.5)  # Returns 1.0
leg(1, 0.5)  # Returns 0.5
leg(2, 0.5)  # Returns -0.125
```
"""
function leg(order::Int, x::Number)
    @assert order >= 0 "Order l must be non-negative"
    T = typeof(x)

    if order == 0
        return one(T)
    elseif order == 1
        return x
    elseif order == 2
        return (T(3) * x^2 - T(1)) / T(2)
    elseif order == 3
        return (T(5) * x^3 - T(3) * x) / T(2)
    elseif order == 4
        return (T(35) * x^4 - T(30) * x^2 + T(3)) / T(8)
    end

    # Start the recurrence with P₀ and P₁.
    previous_legendre = one(T)
    current_legendre = x

    for degree in 1:(order-1)
        next_legendre = ((T(2*degree+1)) * x * current_legendre - T(degree) * previous_legendre) / T(degree+1)
        previous_legendre, current_legendre = current_legendre, next_legendre
    end

    return current_legendre
end

"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Create 2D grid coordinates from 1D coordinate vectors.
Equivalent to MATLAB's `meshgrid` function.

# Arguments
- `x::AbstractVector`: 1D vector of x-coordinates.
- `y::AbstractVector`: 1D vector of y-coordinates.

# Returns
- Tuple `(X, Y)` where:
  - `X` is a matrix of size `(length(y), length(x))` where each row is `x`.
  - `Y` is a matrix of size `(length(y), length(x))` where each column is `y`.

# Examples
```julia
x = 1:3
y = 4:5
X, Y = meshgrid(x, y)
# X = [1 2 3;
#      1 2 3]
# Y = [4 4 4;
#      5 5 5]
```
"""
function meshgrid(x::AbstractVector, y::AbstractVector)
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(reshape(y, :, 1), 1, length(x))
    return X, Y
end

"""
    cart2pol(x, y)

Convert Cartesian coordinates (x, y) to polar coordinates (radius, angle).
Works element-wise for arrays.

# Arguments
- `x`: Scalar or array of x-coordinates.
- `y`: Scalar or array of y-coordinates (must match size of `x` if arrays).

# Returns
- Tuple `(r, θ)` where:
  - `r` is the radius (sqrt(x² + y²)).
  - `θ` is the angle in radians, using `atan(y, x)` for quadrant correctness (range [-π, π]).
"""
function cart2pol(x, y)
    r = sqrt.(x.^2 .+ y.^2)
    θ = atan.(y, x)
    return r, θ
end

"""
    pol2cart(r, θ)

Convert polar coordinates (radius, angle) to Cartesian coordinates (x, y).
Works element-wise for arrays.

# Arguments
- `r`: Scalar or array of radii.
- `θ`: Scalar or array of angles in radians (must match size of `r` if arrays).

# Returns
- Tuple `(x, y)` where:
  - `x = r * cos(θ)`.
  - `y = r * sin(θ)`.
"""
function pol2cart(r, θ)
    x = r .* cos.(θ)
    y = r .* sin.(θ)
    return x, y
end

end # module MathUtils