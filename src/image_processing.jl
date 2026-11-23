module ImageProcessing

export foldQuadrant, unfoldQuadrant, resizeFolded

"""
    foldQuadrant(M, x0, y0; quadrant_filter=[1, 1, 1, 1])

Sum selected quadrants of a detector image or stack around the one-based
center `(x0, y0)` in (column, row) order. The filter follows MATLAB CPBASEX
order: top-right, top-left, bottom-left, bottom-right. Each selected
quadrant includes the center row and column. The result starts at the
center and extends as far as the smallest selected quadrant in each direction.
"""
function foldQuadrant(M::AbstractArray{T,N}, x0::Int, y0::Int;
                      quadrant_filter::Vector{Int}=[1, 1, 1, 1]) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert length(quadrant_filter) == 4 "quadrant_filter must have exactly 4 elements"
    @assert all(in((0, 1)), quadrant_filter) "quadrant_filter must contain only 0s or 1s"
    selected = findall(==(1), quadrant_filter)
    @assert !isempty(selected) "quadrant_filter must select at least one quadrant"

    rows, cols = size(M, 1), size(M, 2)
    @assert 1 <= x0 <= cols "Center x0=$x0 is outside the image columns"
    @assert 1 <= y0 <= rows "Center y0=$y0 is outside the image rows"

    directions = ((-1, 1), (-1, -1), (1, -1), (1, 1))
    lx = minimum(i -> directions[i][2] > 0 ? cols - x0 : x0 - 1, selected)
    ly = minimum(i -> directions[i][1] > 0 ? rows - y0 : y0 - 1, selected)
    trailing = N == 3 ? (Colon(),) : ()
    out = zeros(T, (ly + 1, lx + 1, size(M)[3:end]...))
    for i in selected
        dy, dx = directions[i]
        ys = range(y0; step=dy, length=ly + 1)
        xs = range(x0; step=dx, length=lx + 1)
        out .+= view(M, ys, xs, trailing...)
    end
    return out
end

"""
    unfoldQuadrant(M; center_offset=true, quadrant_count=4)

Mirror a quadrant or stack in both directions and divide by the number of
quadrants summed in `M`. With `center_offset=true`, the first row and column
become the center of an odd-sized output. Set it to `false` to duplicate them.
"""
function unfoldQuadrant(M::AbstractArray{T,N}; center_offset::Bool=true,
                        quadrant_count::Int=4) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert 1 <= quadrant_count <= 4 "quadrant_count must be between 1 and 4"
    first_mirrored = 1 + Int(center_offset)
    rows, cols = size(M, 1), size(M, 2)
    @assert rows >= first_mirrored && cols >= first_mirrored "Input quadrant is too small for center_offset=$center_offset"
    ys = [rows:-1:first_mirrored; 1:rows]
    xs = [cols:-1:first_mirrored; 1:cols]
    trailing = N == 3 ? (Colon(),) : ()
    return M[ys, xs, trailing...] ./ quadrant_count
end

"""
    resizeFolded(data, rmax)

Resize a folded quadrant or stack to `rmax × rmax` by cropping or
zero-padding each spatial dimension. The center stays at index `[1, 1]`.
"""
function resizeFolded(data::AbstractArray{T,N}, rmax::Int) where {T,N}
    @assert N in (2, 3) "Input must be a 2D or 3D array"
    @assert rmax > 0 "Target size rmax must be positive"
    rows = 1:min(size(data, 1), rmax)
    cols = 1:min(size(data, 2), rmax)
    trailing = N == 3 ? (Colon(),) : ()
    out = zeros(T, (rmax, rmax, size(data)[3:end]...))
    out[rows, cols, trailing...] .= view(data, rows, cols, trailing...)
    return out
end

end # module ImageProcessing
