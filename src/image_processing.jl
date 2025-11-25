"""
    ImageProcessing

Provides functions for manipulating image quadrants, specifically for folding multiple
quadrants onto one, unfolding a single quadrant into a symmetric image, and resizing
folded quadrant images. These operations are often used in analyzing data with
cylindrical or four-fold symmetry, such as photoelectron spectroscopy images.
"""
module ImageProcessing

export foldQuadrant, unfoldQuadrant, resizeFolded

"""
    foldQuadrant(M::AbstractArray{T,N}, x0::Int, y0::Int; 
                 quadrant_filter::Vector{Int}=[1,1,1,1]) where {T, N} -> AbstractArray{T,N}

Folds specified image quadrants onto the first quadrant relative to a center point (x0, y0).

This function takes a 2D or 3D image array `M` and adds the pixel values from the selected
quadrants (defined by `quadrant_filter`) into a new array representing the first quadrant region.
The size of the output quadrant is determined by the minimum available space from the center
to the edges across the selected quadrants.

# Arguments
- `M::AbstractArray{T,N}`: Input image array (N=2 for 2D, N=3 for 3D).
- `x0::Int`: x-coordinate (column index) of the folding center.
- `y0::Int`: y-coordinate (row index) of the folding center.
- `quadrant_filter::Vector{Int}`: Optional. A 4-element vector containing 0s or 1s.
  Indicates which quadrants to include in the sum. Default is `[1,1,1,1]` (all quadrants).
  The order corresponds to standard Cartesian quadrants relative to `(x0, y0)` but using
  image coordinates (row, col):
    - `[1]`: Quadrant II (Top-Left: y <= y0, x <= x0)
    - `[2]`: Quadrant I (Top-Right: y <= y0, x > x0)  <- This is different from standard math convention
    - `[3]`: Quadrant IV (Bottom-Right: y > y0, x > x0) <- This is different from standard math convention
    - `[4]`: Quadrant III (Bottom-Left: y > y0, x <= x0)

# Returns
- `AbstractArray{T,N}`: The folded image quadrant. Its size is `(ly+1, lx+1)` for 2D or
  `(ly+1, lx+1, depth)` for 3D, where `lx` and `ly` are the calculated minimum horizontal
  and vertical extents from the center across the included quadrants.

# Details
- The calculation of `lx` and `ly` ensures the output quadrant fits within the smallest
  bounding box available across all included input quadrants relative to the center.
  The `10^10` effectively acts as infinity for excluded quadrants, ensuring they don't
  restrict the size.
- Assumes standard Julia array indexing `M[row, column]`, equivalent to `M[y, x]`.
"""
function foldQuadrant(M::AbstractArray{T,N}, x0::Int, y0::Int; 
                      quadrant_filter::Vector{Int}=[1,1,1,1]) where {T, N}
    # Validate inputs
    @assert N == 2 || N == 3 "Input must be a 2D or 3D array, got $(N)D"
    @assert length(quadrant_filter) == 4 "quadrant_filter must have exactly 4 elements"
    @assert all(val -> val == 0 || val == 1, quadrant_filter) "quadrant_filter must contain only 0s or 1s"

    # Get dimensions (rows = y, cols = x)
    rows, cols = size(M, 1), size(M, 2)

    # Check center coordinates are valid (1-based indexing)
    @assert 1 <= x0 <= cols "Center x0=$x0 is outside image column bounds [1, $(cols)]"
    @assert 1 <= y0 <= rows "Center y0=$y0 is outside image row bounds [1, $(rows)]"

    depth = (N == 3) ? size(M, 3) : 1

    # Calculate the maximum extent (lx, ly) for the output quadrant.
    # This depends on the distance from the center (x0, y0) to the edges
    # in each direction, considering only the quadrants included by the filter.
    # A large number (1e10) is used to effectively ignore distances for filtered-out quadrants.
    large_num = 1e10 # Represents infinity for filtering purposes
    lx = floor(Int, min(
        # Max extent towards left (cols - x0) vs right (x0 - 1), considering filters
        max(cols - x0, min((1 - quadrant_filter[1]) * large_num, (1 - quadrant_filter[4]) * large_num)), # Max extent in +x direction (Quad I, IV)
        max(x0 - 1,    min((1 - quadrant_filter[2]) * large_num, (1 - quadrant_filter[3]) * large_num))  # Max extent in -x direction (Quad II, III)
    ))
    ly = floor(Int, min(
        # Max extent downwards (rows - y0) vs upwards (y0 - 1), considering filters
        max(rows - y0, min((1 - quadrant_filter[3]) * large_num, (1 - quadrant_filter[4]) * large_num)), # Max extent in +y direction (Quad III, IV)
        max(y0 - 1,    min((1 - quadrant_filter[1]) * large_num, (1 - quadrant_filter[2]) * large_num))  # Max extent in -y direction (Quad I, II)
    ))

    # Initialize output array (size ly+1, lx+1 because ranges include the center)
    if N == 3
        M_out = zeros(T, ly + 1, lx + 1, depth)
    else
        M_out = zeros(T, ly + 1, lx + 1)
    end

    # Define sign multipliers for accessing each quadrant relative to (x0, y0)
    # Format: [sign_y, sign_x] for [Quad II, Quad I, Quad IV, Quad III]
    # Note the quadrant order [II, I, IV, III] matches the filter index [1, 2, 3, 4]
    quad_signs = [-1  1; # Q1 (y0-ly : y0, x0 : x0+lx) -> Top-Right
                  -1 -1; # Q2 (y0-ly : y0, x0-lx : x0) -> Top-Left
                   1 -1; # Q3 (y0 : y0+ly, x0-lx : x0) -> Bottom-Left
                   1  1] # Q4 (y0 : y0+ly, x0 : x0+lx) -> Bottom-Right

    # Process each quadrant
    for i in 1:4
        if quadrant_filter[i] == 1
            sign_y, sign_x = quad_signs[i, 1], quad_signs[i, 2]

            # Define ranges for reading from the input image M
            # These ranges select the appropriate quadrant relative to (x0, y0)
            # The step direction depends on the signs.
            y_range = range(y0, step=sign_y, length=ly+1)
            x_range = range(x0, step=sign_x, length=lx+1)

            # Extract the quadrant data
            # Note: Indexing M[y_range, x_range] correctly handles reversed ranges
            quadrant_data = if N == 3
                M[y_range, x_range, :]
            else
                M[y_range, x_range]
            end
            M_out .+= quadrant_data
        end
    end

    return M_out
end

"""
    unfoldQuadrant(M::AbstractArray{T,N}; center_offset::Bool=true) where {T, N} -> AbstractArray{T,N}

Unfolds a quadrant (assumed to be the first quadrant) into a full image with fourfold symmetry.

This function takes an array `M` representing a single quadrant (typically the result of
`foldQuadrant` or a similar operation) and mirrors it horizontally and vertically
to create a complete, symmetric image. It's the inverse operation of folding all four
quadrants into one.

# Arguments
- `M::AbstractArray{T,N}`: Input folded quadrant array (N=2 for 2D, N=3 for 3D).
- `center_offset::Bool`: Optional. Determines how the center row/column is handled during mirroring.
    - `true` (default): Assumes the input quadrant shares the center row/column with adjacent
      quadrants (typical for images with an even total size after unfolding). The first row/column
      of `M` is treated as the center and is *not* duplicated in the mirroring.
    - `false`: Assumes the input quadrant represents a unique center (typical for images with
      an odd total size after unfolding). The first row/column of `M` is duplicated during mirroring.

# Returns
- `AbstractArray{T,N}`: The full unfolded image with fourfold symmetry.
  The output dimensions will be approximately `2*size(M,1) x 2*size(M,2)` (minus 1 row/col if `center_offset=true`).
  The pixel values are normalized by dividing by the number of contributing quadrants (typically 4 if
  `M` was produced by folding all four quadrants).

# Details
- The function constructs the full image by concatenating the original quadrant (`M`)
  with its flipped versions.
- Division by 4 assumes `M` was created by summing 4 quadrants; adjust if needed.
"""
function unfoldQuadrant(M::AbstractArray{T,N}; center_offset::Bool=true) where {T, N}
    # Validate input
    @assert N == 2 || N == 3 "Input must be a 2D or 3D array, got $(N)D"

    # Determine the offset for slicing based on center_offset
    # offset=1 means skip the first row/col when flipping (center is shared)
    # offset=0 means include the first row/col when flipping (center is unique)
    offset = center_offset ? 1 : 0
    rows, cols = size(M, 1), size(M, 2)

    # Check if there's enough data to mirror based on the offset
    @assert rows > offset "Input quadrant height ($(rows)) too small for center_offset=$(center_offset)"
    @assert cols > offset "Input quadrant width ($(cols)) too small for center_offset=$(center_offset)"

    # Select the parts to be flipped (excluding the center row/col if offset=1)
    M_flip_rows = 1+offset:rows
    M_flip_cols = 1+offset:cols

    # Mirroring logic (handles 2D and 3D)
    if N == 3
        # Top-Left Quadrant (flipped vertically and horizontally)
        TL = M[M_flip_rows, M_flip_cols, :] |> arr -> reverse(arr, dims=1) |> arr -> reverse(arr, dims=2)
        # Top-Right Quadrant (flipped vertically)
        TR = M[M_flip_rows, :, :] |> arr -> reverse(arr, dims=1)
        # Bottom-Left Quadrant (flipped horizontally)
        BL = M[:, M_flip_cols, :] |> arr -> reverse(arr, dims=2)
        # Bottom-Right Quadrant (original)
        BR = M

        # Concatenate quadrants
        M_out = cat(
            cat(TL, TR; dims=2), # Top row
            cat(BL, BR; dims=2);  # Bottom row
            dims=1
        )
    else # N == 2
        # Top-Left Quadrant (flipped vertically and horizontally)
        TL = M[M_flip_rows, M_flip_cols] |> arr -> reverse(arr, dims=1) |> arr -> reverse(arr, dims=2)
        # Top-Right Quadrant (flipped vertically)
        TR = M[M_flip_rows, :] |> arr -> reverse(arr, dims=1)
        # Bottom-Left Quadrant (flipped horizontally)
        BL = M[:, M_flip_cols] |> arr -> reverse(arr, dims=2)
        # Bottom-Right Quadrant (original)
        BR = M

        # Concatenate quadrants using matrix notation
        M_out = [ TL TR ; 
                  BL BR ]
    end

    # Normalize: Assuming M was created by summing 4 quadrants.
    # If M is just a single quadrant measurement, this division might not be desired.
    return M_out ./ 4
end

"""
    resizeFolded(data::AbstractArray{T,N}, rmax::Int) where {T, N} -> AbstractArray{T,N}

Resizes a square, folded quadrant image to a new square size `(rmax, rmax)`.

The input `data` is assumed to be a square array representing the first quadrant
(e.g., from `foldQuadrant`). The resizing is done by either cropping or padding
with zeros, anchored at the top-left corner (index [1, 1]).

# Arguments
- `data::AbstractArray{T,N}`: Input square folded quadrant (N=2 for 2D, N=3 for 3D).
- `rmax::Int`: The desired side length for the output square array.

# Returns
- `AbstractArray{T,N}`: The resized square array of size `(rmax, rmax)` or `(rmax, rmax, depth)`.
  - If `size(data, 1) > rmax`, the output is the top-left `rmax x rmax` portion of `data`.
  - If `size(data, 1) < rmax`, the output is `data` padded with zeros at the bottom and right
    to reach the size `rmax x rmax`.

# Errors
- Throws an `ArgumentError` if the input `data` is not square (rows != cols).
"""
function resizeFolded(data::AbstractArray{T,N}, rmax::Int) where {T, N}
    # Validate input
    @assert rmax > 0 "Target size rmax must be positive, got $rmax"
    @assert N == 2 || N == 3 "Input must be 2D or 3D array, got $(N)D"

    rows = size(data, 1)
    cols = size(data, 2)
    @assert rows == cols "Input data must be square (rows == cols), got ($rows, $cols)"
    current_size = rows
    depth = (N == 3) ? size(data, 3) : 1

    # Determine crop/pad dimensions
    if current_size == rmax
        return copy(data) # No change needed, return a copy
    elseif current_size > rmax
        # Crop: Select the top-left rmax x rmax part
        if N == 3
            return data[1:rmax, 1:rmax, :]
        else
            return data[1:rmax, 1:rmax]
        end
    else # current_size < rmax
        # Pad: Create a larger array and place data in the top-left
        pad_rows = rmax - current_size
        pad_cols = rmax - current_size

        if N == 3
            # Create the padded output array initialized with zeros
            out = zeros(T, rmax, rmax, depth)
            # Place the original data in the top-left corner
            out[1:current_size, 1:current_size, :] = data
            return out
        else # N == 2
            # Create the padded output array initialized with zeros
            out = zeros(T, rmax, rmax)
            # Place the original data in the top-left corner
            out[1:current_size, 1:current_size] = data
            return out
        end
    end
end

end # module ImageProcessing 