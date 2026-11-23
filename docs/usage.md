# Run an inversion

Start Julia with `julia --project=. --threads=auto` from the repository root.
Use a square `Float64` or `Float32` matrix of finite intensities whose rows and
columns follow the detector image. One `pbasex` call processes one image.

```julia
using pBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
result = pbasex(image; center=(512, 512), make_images=false)
```

Set `make_images=true` to include the fitted detector projection and the
Abel-inverted central slice. See [basis selection](bases.md#select-a-basis)
to reuse an existing Julia or MATLAB basis or to control where a new one
is saved.

## Set the image center and quadrants

`center` is a one-based integer tuple `(column, row)`. By default, Julia uses
`(div(size(image, 2), 2) + 1, div(size(image, 1), 2) + 1)`. Supply the measured
center when it differs from this default. Convert a NumPy zero-based
`[row, column]` pair by reversing its order and adding one to both coordinates.
For example, `[473, 476]` becomes `(477, 474)`.

`quadrant_filter` takes a four-element integer vector of zeros and ones in
CPBASEX order.

| Index | Quadrant |
| --- | --- |
| 1 | Top-right |
| 2 | Top-left |
| 3 | Bottom-left |
| 4 | Bottom-right |

Use `[1, 1, 1, 1]` for all quadrants or `[1, 0, 0, 0]` for only the top-right.
The fit scales the selected sum by four divided by the number selected. For
a symmetric image, choosing fewer quadrants therefore preserves the
four-quadrant intensity scale.

The folded quadrant starts at the center and extends as far as the selected
quadrants have pixels in common. Julia crops or zero-pads each dimension to
the basis size. With all quadrants selected, the last radius with complete
detector coverage is the shortest distance from the center to any image edge.
For the 1024 × 1024 example centered at `(477, 474)`, this radius is 473 px.
Bins at larger radii depend on a fit to incomplete detector coverage.

## Set the energy scale and fit

| Keyword | Default | Meaning |
| --- | --- | --- |
| `alpha` | `1.0` | Positive calibration factor in `E = alpha * r^2`. With the default, the axis is in pixel². |
| `reg` | `0.0` | Nonnegative Tikhonov regularization. With positive `reg`, the solver uses `s / (s^2 + reg)` for singular value `s`. At zero it takes a pseudoinverse with a numerical cutoff. |
| `make_images` | `true` | Return the fitted projection and inverted slice. Set to `false` to skip generating or reading the central-slice matrix. |
| `weights` | `nothing` | Finite pixel multipliers with the same size and element type as the input image. |

`IE` is the raw fitted one-dimensional energy distribution. The radial to
energy conversion multiplies the isotropic radial profile by `r / (2 * alpha)`.
Changing `alpha` therefore changes both the energy axis and intensity per unit
energy. If you integrate the spectrum, use the spacing in `E` because the bins
are unevenly spaced on that axis.

### Weighted fits

Julia folds and resizes the weight image separately, then multiplies both the
folded detector image and the corresponding rows of the forward matrix by
the processed weights. For a desired inverse-variance weight in the
least-squares objective, use an inverse standard deviation as the multiplier.
Keep the original weight image, center, and quadrant selection together when
reusing a named weighted basis. Reconstructed projection pixels with zero
folded weight are `NaN`.

## Read and save the results

| Field | Shape and meaning |
| --- | --- |
| `E` | Vector of `length(x)` energy samples on the basis grid. |
| `IE` | Vector of `length(x)` intensity values. |
| `betas` | Matrix with one row per energy and one column per nonzero Legendre order, in basis order. At zero `IE`, entries are `NaN`. |
| `coefficients` | Matrix with `length(k)` radial centers as rows and `length(l)` Legendre orders as columns. |
| `reconstructed_image` | Fitted detector projection, or `nothing` with `make_images=false`. |
| `inverted_image` | Abel-inverted central slice, or `nothing` with `make_images=false`. |
| `residuals` | Input minus reconstructed projection when both grids have the same size and center. Otherwise `nothing`. |

An image reconstructed from a basis of length `nx` has side length
`2 * nx - 1` and center `(nx, nx)`. This makes the output grid independent of
the input image size and measured center.

Save the spectrum and anisotropies from the example above using HDF5.

```julia
using HDF5

h5open("results.h5", "w") do file
    write(file, "E", result.E)
    write(file, "IE", result.IE)
    write(file, "betas", result.betas)
end
```

For a sequence of images, loop over matrices and pass the same named basis
when they share a grid and weighting. The same basis file can be used with
different values of `alpha` and `reg`.
