# Run an inversion

Start Julia with `julia --project=. --threads=auto` from the repository root.
Use a `Float64` or `Float32` matrix of finite intensities whose rows and
columns follow the detector image. One `pbasex` call processes one image.
The default quadrant geometry requires a square image. Half-plane
geometry also accepts rectangular images.

```julia
using PBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
basis_path = joinpath(mkpath("data"), "gaussian_r512_k128_l4.h5")
result = pbasex(image; basis=basis_path, center=(512, 512), make_images=false)
```

The required `basis` keyword takes a file path. Relative paths start from
the Julia working directory. PBASEX reads an existing file or creates one
at a new path whose parent directory exists. Set `make_images=true` to
include the fitted detector projection and Abel-inverted central slice.
[Basis selection](bases.md#select-a-basis)

| Geometry | Detector symmetry used | Default orders for a new basis |
| --- | --- | --- |
| `:quadrant` | Left–right and up–down reflection | `0, 2, 4` |
| `:half_plane` | Left–right reflection about the vertical physical axis | `0, 1, 2, 3, 4` |

## Set the image center and quadrants

`center` is a one-based integer tuple `(column, row)`. By default, Julia uses
`(div(size(image, 2), 2) + 1, div(size(image, 1), 2) + 1)`. Supply the measured
center when it differs from this default. Convert a NumPy zero-based
`[row, column]` pair by reversing its order and adding one to both coordinates.
For example, `[473, 476]` becomes `(477, 474)`.

### Align the polarization axis

The angular basis uses the vertical detector coordinate to calculate
`cos(theta)`. Supply an image with its polarization axis vertical. For a
horizontal axis, rotate the image 90° counterclockwise with `rotl90` and
transform its measured center using the original image width.

```julia
# For a horizontal-polarization image with NumPy center [473, 476].
column, row = 477, 474
aligned_image = rotl90(image)
aligned_center = (row, size(image, 2) + 1 - column)
result = pbasex(aligned_image; basis=basis_path, center=aligned_center, make_images=false)
```

For a 1024 × 1024 image, this gives `(474, 548)` in one-based
`(column, row)` coordinates. A transpose gives `(474, 477)` because it
swaps axes without reversing the new row axis. With all four quadrants
and even Legendre orders, the transposed and rotated images give
equivalent fits. With odd orders, the choice reverses the positive
direction in this example and therefore the odd beta signs. Quadrant
selection refers to the transformed image.

### Fit forward–backward asymmetry

Use `geometry=:half_plane` for a cylindrically symmetric image whose
vertical hemispheres can differ. Align the physical axis vertically and
orient its positive direction toward the top of the image. Pass a separate
basis filename for this geometry.

```julia
signed_basis = joinpath(mkpath("data"), "gaussian-half-plane-r512-l4.h5")
result = pbasex(aligned_image; basis=signed_basis, center=aligned_center,
                geometry=:half_plane, l=0:4, make_images=false)
result.legendre_orders[2:end]  # The order of the beta columns
```

Half-plane mode averages the selected left and right pixels at each signed
axial coordinate. The axis column contributes once. Reconstructions mirror
only the columns. Flipping the image top to bottom and updating its center
reverses odd beta signs. The angle-integrated spectrum still comes from
order zero. Both geometries fit in detector-pixel intensity units.

The selected quadrants must include at least one from each hemisphere.
Where only one selected side has a pixel, the fit uses that measurement.
The half-plane fit excludes folded positions with no measured pixel.
Reuse the same detector coverage when reopening a named half-plane basis.
[Theory and analytic example](forward-backward-asymmetry.md)

### Select quadrants

`quadrant_filter` takes a four-element integer vector of zeros and ones in
CPBASEX order.

| Index | Quadrant |
| --- | --- |
| 1 | Top-right |
| 2 | Top-left |
| 3 | Bottom-left |
| 4 | Bottom-right |

Use `[1, 1, 1, 1]` for all quadrants. In quadrant mode, `[1, 0, 0, 0]`
selects only the top-right. The quadrant fit averages the selected quadrants.
It infers the full distribution using the assumed reflection symmetries.
In half-plane mode, `[1, 0, 0, 1]`
selects the right side in both hemispheres.

In quadrant mode, the folded quadrant starts at the center and extends as far as the selected
quadrants have pixels in common. Julia crops or zero-pads each dimension to
the basis size. With all quadrants selected, the last radius with complete
detector coverage is the shortest distance from the center to any image edge.
For the aligned 1024 × 1024 example centered at `(474, 548)`, this radius
is 473 px.
Bins at larger radii depend on a fit to incomplete detector coverage.

## Set the energy scale and fit

| Keyword | Default | Meaning |
| --- | --- | --- |
| `alpha` | `1.0` | Positive calibration factor in `E = alpha * r^2`. With the default, the axis is in pixel². |
| `reg` | `0.0` | Nonnegative Tikhonov regularization. With positive `reg`, the solver uses `s / (s^2 + reg)` for singular value `s`. At zero it takes a pseudoinverse with a numerical cutoff. |
| `make_images` | `true` | Return the fitted projection and inverted slice. Set to `false` to skip generating or reading the central-slice matrix. |

### Physical units and peak yields

For a histogram of detected electrons on a unit-pixel grid, both geometries
return the following quantities. Any prior rescaling of the input intensity
carries through to the spectrum and both images.

- `IE` samples the angle-integrated yield density $dN/dE$. Its units are
  electrons/eV when $E$ is calibrated in eV, or electrons/px² at `alpha=1`.
  The fitted electron yield in a peak is $\int_{E_a}^{E_b} I(E)\,dE$.
  Integrate with the actual spacing of `result.E`, interpolating at the
  window edges. The factor $r/(2\alpha)$ in `IE` already accounts for the
  change of variable $E=\alpha r^2$. The least-squares fit does not impose
  exact count conservation.
- `reconstructed_image` gives the fitted detector projection in electrons
  per detector pixel. Each pixel samples the line-of-sight integral of the
  fitted 3D density. `residuals` uses the same units.
- `inverted_image` samples that 3D density in the central plane, in
  electrons/px³ for velocity coordinates expressed in detector pixels.
  A count in a finite slab requires integrating the density over the slab
  volume. Summing the 2D slice alone omits its depth.

For raw detector counts, `pbasex_counts` performs one inversion and propagates
counting errors through the fit. Choose the energy intervals for the observable
you want to report.

| Observable | What to pass | Where to read the values and errors |
| --- | --- | --- |
| [Spectrum and betas at each native energy](counting-errors.md#spectrum-and-betas-on-the-native-grid) | Omit `windows` | `analysis.result.IE` and `analysis.spectrum_se`, or `analysis.result.betas` and `analysis.beta_se` |
| [Integrated yield and betas for selected peaks](counting-errors.md#integrated-yield-and-betas-for-selected-peaks) | One `(low, high)` pair per interval in `windows` | `analysis.windows[i]` contains `yield`, `yield_se`, `betas`, and `beta_se` |
| [Bin-averaged spectrum and betas on a uniform grid](counting-errors.md#yield-preserving-spectrum-and-betas-on-a-uniform-energy-grid) | Adjacent, equal-width intervals in `windows` | Read each bin's `spectrum` and `spectrum_se`. Use its betas and beta errors directly. |

## Read and save the results

For a corrected event checkpoint with shot IDs, use `pbasex_shots` to
estimate uncertainty from the spread across whole shots. It accepts an
existing basis and any collection of energy windows.

Both `pbasex_counts` and `pbasex_shots` take `center=(column, row)` in
one-based bins of the original histogram. For horizontal polarization,
`rotation=90` turns the detector's right edge toward the top of the fit
image. Use `rotation=-90` to turn its left edge toward the top instead.
This choice sets the positive axial direction for odd Legendre orders.
The default is zero degrees. Multiples of 90° preserve the pixel counts
exactly. Both functions transform the center along with the data.

[Shot-resolved example and output fields](counting-errors.md#shot-resolved-uncertainties-from-a-corrected-checkpoint)

For a uniform spectrum, the [side-by-side Julia example](shot-uncertainties.md#a-uniform-spectrum-from-either-input)
uses the same window definitions and plotting fields for both inputs.
The [step-by-step mathematical derivation](uncertainty-derivation.md) explains
why the fitted yields agree and how each function calculates an error bar.

The following fields belong to the ordinary inversion result. Access them
directly after `pbasex`, or through `analysis.result` after `pbasex_counts`.

| Field | Shape and meaning |
| --- | --- |
| `E` | Vector of `length(x)` energy samples on the basis grid. |
| `IE` | Vector of `length(x)` samples of $dN/dE$. |
| `betas` | Matrix with one row per energy and one column per nonzero Legendre order, in basis order. Each value is the corresponding angular component divided by `IE`. At zero `IE`, entries are `NaN`. |
| `legendre_orders` | Legendre orders in coefficient-column order, including zero. `legendre_orders[2:end]` labels the beta columns. |
| `coefficients` | Matrix with `length(k)` radial centers as rows and `length(l)` Legendre orders as columns. |
| `reconstructed_image` | Fitted detector projection, or `nothing` with `make_images=false`. |
| `inverted_image` | Abel-inverted central slice, or `nothing` with `make_images=false`. |
| `residuals` | Input minus reconstructed projection when both grids have the same size and center. Otherwise `nothing`. |

An image reconstructed from a basis of length `nx` has side length
`2 * nx - 1` and center `(nx, nx)`. This makes the output grid independent of
the input image size and measured center.

Beta values can be sensitive to small fit differences where the intensity
is weak. For a basis comparison, report the intensity range used to
compare the beta curves. [Worked comparison](matlab-migration.md#circularized-image-with-aligned-polarization)

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
when they share a grid and detector coverage. The same basis file can be used with
different values of `alpha` and `reg`.
