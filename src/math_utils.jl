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
function leg(l::Int, x::AbstractArray)
    @assert l >= 0 "Order l must be non-negative"
    T = eltype(x) # Get the element type for output arrays

    # Use explicit formulas for low orders (more efficient and stable)
    if l == 0
        return ones(T, size(x)) # P_0(x) = 1
    elseif l == 1
        return T.(copy(x)) # P_1(x) = x (ensure type T)
    elseif l == 2
        return @. (3 * x^2 - 1) / 2 # P_2(x) = (3x^2 - 1)/2
    elseif l == 3
        return @. (5 * x^3 - 3 * x) / 2 # P_3(x) = (5x^3 - 3x)/2
    elseif l == 4
        return @. (35 * x^4 - 30 * x^2 + 3) / 8 # P_4(x) = (35x^4 - 30x^2 + 3)/8
    end
    
    # Pre-allocate arrays for the recurrence relation
    p_nm2 = ones(T, size(x))
    p_nm1 = T.(copy(x))
    p_n = similar(x, T)
    
    # Apply the recurrence relation: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
    for n in 1:(l-1)
        # Use in-place broadcasting for better performance
        p_n .= ((2*n+1) .* x .* p_nm1 .- n .* p_nm2) ./ (n+1)
        
        # Swap references for next iteration - no allocations
        p_nm2, p_nm1, p_n = p_nm1, p_n, p_nm2
    end
    
    # With the swapping pattern, the result is in p_nm1 after the loop
    return p_nm1
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
function leg(l::Int, x::Number)
    @assert l >= 0 "Order l must be non-negative"
    T = typeof(x) # Preserve the input number type

    # Explicit formulas for low orders
    if l == 0
        return one(T) # P_0(x) = 1
    elseif l == 1
        return x # P_1(x) = x
    elseif l == 2
        return (T(3) * x^2 - T(1)) / T(2) # P_2(x) = (3x^2 - 1)/2
    elseif l == 3
        return (T(5) * x^3 - T(3) * x) / T(2) # P_3(x) = (5x^3 - 3x)/2
    elseif l == 4
        return (T(35) * x^4 - T(30) * x^2 + T(3)) / T(8) # P_4(x) = (35x^4 - 30x^2 + 3)/8
    end

    # For l > 4, use the recurrence relation iteratively
    # Initialize P_0 and P_1
    p_nm2 = one(T) # P_{n-2}, starts with P_0
    p_nm1 = x      # P_{n-1}, starts with P_1

    # Iterate from n=1 up to n=l-1
    for n in 1:(l-1)
        # Calculate P_{n+1} using P_n (p_nm1) and P_{n-1} (p_nm2)
        p_n = ((T(2*n+1)) * x * p_nm1 - T(n) * p_nm2) / T(n+1)
        # Update values for the next iteration: P_{n-1} becomes P_n, P_{n-2} becomes P_{n-1}
        p_nm2, p_nm1 = p_nm1, p_n
    end

    # After the loop finishes (n=l-1), p_nm1 holds P_l(x)
    return p_nm1
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
    # Reshape x into a 1-row matrix and repeat it length(y) times vertically.
    X = repeat(reshape(x, 1, :), length(y), 1)
    # Reshape y into a 1-column matrix and repeat it length(x) times horizontally.
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
    # Calculate radius r = sqrt(x^2 + y^2)
    r = sqrt.(x.^2 .+ y.^2)
    # Calculate angle theta using atan(y, x) for correct quadrant
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
    # Calculate x = r * cos(theta)
    x = r .* cos.(θ)
    # Calculate y = r * sin(theta)
    y = r .* sin.(θ)
    return x, y
end

end # module MathUtils