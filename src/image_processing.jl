module ImageProcessing

export fold_quadrants, unfold_quadrant, resize_folded_quadrant

"""
    fold_quadrants(detector_image, center_column, center_row; quadrant_filter=[1, 1, 1, 1])

Sum selected quadrants of a detector image or stack around the one-based
center `(center_column, center_row)` in (column, row) order. The filter follows MATLAB CPBASEX
order: top-right, top-left, bottom-left, bottom-right. Each selected
quadrant includes the center row and column. The result starts at the
center and extends as far as the smallest selected quadrant in each direction.
"""
function fold_quadrants(detector_image::AbstractArray{T,N}, center_column::Int, center_row::Int;
                        quadrant_filter::Vector{Int}=[1, 1, 1, 1]) where {T,N}
    return quadrant_sum(detector_image, center_column, center_row, quadrant_filter, false)
end

function quadrant_sum(detector_image::AbstractArray{T,N}, center_column, center_row,
                      quadrant_filter, poisson_variance) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert length(quadrant_filter) == 4 "quadrant_filter must have exactly 4 elements"
    @assert all(in((0, 1)), quadrant_filter) "quadrant_filter must contain only 0s or 1s"
    selected_quadrants = findall(==(1), quadrant_filter)
    @assert !isempty(selected_quadrants) "quadrant_filter must select at least one quadrant"

    image_rows, image_columns = size(detector_image, 1), size(detector_image, 2)
    @assert 1 <= center_column <= image_columns "Center x0=$center_column is outside the image columns"
    @assert 1 <= center_row <= image_rows "Center y0=$center_row is outside the image rows"

    # Steps are ordered as (row, column) to match array indexing.
    quadrant_directions = ((-1, 1), (-1, -1), (1, -1), (1, 1))
    column_extent = minimum(quadrant -> quadrant_directions[quadrant][2] > 0 ? image_columns - center_column : center_column - 1, selected_quadrants)
    row_extent = minimum(quadrant -> quadrant_directions[quadrant][1] > 0 ? image_rows - center_row : center_row - 1, selected_quadrants)
    frame_indices = N == 3 ? (Colon(),) : ()
    folded_projection = zeros(T, (row_extent + 1, column_extent + 1, size(detector_image)[3:end]...))
    for quadrant in selected_quadrants
        row_step, column_step = quadrant_directions[quadrant]
        row_indices = range(center_row; step=row_step, length=row_extent + 1)
        column_indices = range(center_column; step=column_step, length=column_extent + 1)
        pixels = view(detector_image, row_indices, column_indices, frame_indices...)
        folded_projection .+= pixels
        if poisson_variance
            # Each occurrence of the same count contributes its multiplicity
            # times that count to Var(sum), giving multiplicity² in total.
            axial_partner = (2, 1, 4, 3)[quadrant]
            equatorial_partner = (4, 3, 2, 1)[quadrant]
            if quadrant_filter[axial_partner] == 1
                @views folded_projection[2:end, 1] .+= pixels[2:end, 1]
            end
            if quadrant_filter[equatorial_partner] == 1
                @views folded_projection[1, 2:end] .+= pixels[1, 2:end]
            end
            folded_projection[1, 1] += (length(selected_quadrants) - 1) * pixels[1, 1]
        end
    end
    return folded_projection
end

"""
    unfold_quadrant(folded_projection; center_offset=true, quadrant_count=4)

Mirror a quadrant or stack in both directions and divide by the number of
quadrants summed in `folded_projection`. With `center_offset=true`, the first row and column
become the center of an odd-sized output. Set it to `false` to duplicate them.
"""
function unfold_quadrant(folded_projection::AbstractArray{T,N}; center_offset::Bool=true,
                         quadrant_count::Int=4) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert 1 <= quadrant_count <= 4 "quadrant_count must be between 1 and 4"
    first_reflected_index = 1 + Int(center_offset)
    folded_rows, folded_columns = size(folded_projection, 1), size(folded_projection, 2)
    @assert folded_rows >= first_reflected_index && folded_columns >= first_reflected_index "Input quadrant is too small for center_offset=$center_offset"
    row_indices = [folded_rows:-1:first_reflected_index; 1:folded_rows]
    column_indices = [folded_columns:-1:first_reflected_index; 1:folded_columns]
    frame_indices = N == 3 ? (Colon(),) : ()
    return folded_projection[row_indices, column_indices, frame_indices...] ./ quadrant_count
end

"""
    resize_folded_quadrant(folded_projection, detector_side)

Resize a folded quadrant or stack to `detector_side × detector_side` by cropping or
zero-padding each spatial dimension. The center stays at index `[1, 1]`.
"""
function resize_folded_quadrant(folded_projection::AbstractArray{T,N}, detector_side::Int) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert detector_side > 0 "Target size rmax must be positive"
    row_indices = 1:min(size(folded_projection, 1), detector_side)
    column_indices = 1:min(size(folded_projection, 2), detector_side)
    frame_indices = N == 3 ? (Colon(),) : ()
    resized_projection = zeros(T, (detector_side, detector_side, size(folded_projection)[3:end]...))
    resized_projection[row_indices, column_indices, frame_indices...] .= view(folded_projection, row_indices, column_indices, frame_indices...)
    return resized_projection
end

# Keep both axial hemispheres in detector row order. Mirror only the columns.
function fold_half_plane(image::AbstractMatrix{T}, center, detector_side, quadrant_filter, weights;
                         count_variance=nothing) where T
    rows, columns = size(image)
    center_column, center_row = center
    1 <= center_column <= columns && 1 <= center_row <= rows ||
        throw(ArgumentError("Center must lie inside the detector image"))
    length(quadrant_filter) == 4 && all(in((0, 1)), quadrant_filter) ||
        throw(ArgumentError("quadrant_filter must have four entries, each 0 or 1"))
    top_right, top_left, bottom_left, bottom_right = isone.(quadrant_filter)
    (top_left || top_right) && (bottom_left || bottom_right) ||
        throw(ArgumentError("Half-plane fitting requires a selected quadrant in each hemisphere"))

    folded = zeros(T, 2 * detector_side - 1, detector_side)
    multipliers = zeros(T, size(folded))
    for row in axes(folded, 1)
        image_row = center_row + row - detector_side
        1 <= image_row <= rows || continue
        use_left = row < detector_side ? top_left : row > detector_side ? bottom_left : top_left || bottom_left
        use_right = row < detector_side ? top_right : row > detector_side ? bottom_right : top_right || bottom_right
        for column in axes(folded, 2)
            offset = column - 1
            precision_sum = zero(T)
            weighted_signal = zero(T)
            for (direction, selected) in ((1, use_right), (-1, use_left))
                selected || continue
                iszero(offset) && direction == -1 && use_right && continue
                image_column = center_column + direction * offset
                1 <= image_column <= columns || continue
                precision = isnothing(weights) ? one(T) : abs2(weights[image_row, image_column])
                iszero(precision) && continue
                weighted_signal += precision * image[image_row, image_column]
                precision_sum += precision
            end
            if !iszero(precision_sum)
                folded[row, column] = weighted_signal / precision_sum
                if !isnothing(count_variance)
                    # Raw counts, with each distinct selected pixel used once.
                    count_variance[row, column] = folded[row, column] / precision_sum
                end
                # The unweighted mode gives every observed folded pixel equal weight.
                multipliers[row, column] = isnothing(weights) ? one(T) : sqrt(precision_sum)
            end
        end
    end
    return folded, isnothing(weights) && all(isone, multipliers) ? nothing : multipliers
end

function unfold_half_plane(folded)
    columns = size(folded, 2)
    return folded[:, [columns:-1:2; 1:columns]]
end

end # module ImageProcessing
