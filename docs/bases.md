# Build and save a basis

A basis combines radial profiles with Legendre polynomials. The forward
matrix projects each basis function onto the detector. The solver stores
its singular value decomposition and uses the same fit coefficients to
evaluate the radial profiles and central slice.

## Select a basis

| `basis` argument to `pbasex` | Behavior |
| --- | --- |
| Existing HDF5 file | Read the grid and factors from that file. Explicitly supplied settings must agree with the file's metadata. Julia and MATLAB CPBASEX layouts are supported. |
| New HDF5 path | Compute and save the basis at that path. The parent directory must exist. |

Choose a filename that identifies the definition you used. For example,
`gaussian_r512_k256_l10_sigma1p6.h5` identifies a 512-point detector grid,
256 radial centers, even Legendre orders through 10, and a Gaussian width
of 1.6 px. PBASEX reads the definition from the HDF5 file.

```julia
using PBASEX

basis_path = joinpath(mkpath("data"), "gaussian_r512_k256_l10_sigma1p6.h5")
basis_data = BASISData(basis_path; x=0:511, k=0.5:2:510.5, l=0:2:10, params=[1.6])
```

Pass that path to `pbasex(image; basis=basis_path)` after loading an image.
An existing file supplies its basis settings. If you also pass settings
that conflict with the file, PBASEX raises an `ArgumentError`. Use a new
path when you change the basis definition. You can change `alpha` and
`reg` while reusing the same file.

Use `BASISData(; x, k, l)` to build in memory. Both constructors default to
`Float64`. Use `BASISData{Float32}` for single precision.

### Migrate an existing script

Replace an omitted `basis` argument or a basis directory with the path to
one HDF5 file. You can use a previously generated hash-named file by passing
its full path. When loading a MATLAB basis, let the file supply the grid
and matrices. Omit `params` or `radial_basis` if the file cannot verify them.
[MATLAB basis files](matlab-migration.md#reuse-a-matlab-basis-in-julia)

## Set the basis definition

| Keyword | Default for a new basis | Meaning |
| --- | --- | --- |
| `x` | `0:511` | Integer radial grid. Use `0:rmax` for a detector sampled at unit pixel spacing. |
| `k` | `x[2:4:end] .+ 0.5` | Radial centers. |
| `l` | `[0, 2, 4]` | Distinct nonnegative Legendre orders, starting with zero for image inversion. Quadrant fitting uses even orders. |
| `geometry` | `:quadrant` | Use `:half_plane` to retain both signs of the vertical axis and fit even and odd orders. |
| `radial_basis` | `:gaussian` | Built-in family. Use `:cubic_bspline` for compact support. |
| `params` | Depends on the family | First value is the Gaussian width or the B-spline knot spacing, in pixels. |
| `rBF` | Selected family | Custom radial function with signature `rBF(r, k, params)`. |
| `zIP` | Selected family | Custom upper limit with signature `zIP(R, k, params)`, where `R` is the projected radius. |
| `cache_key` | `nothing` | Optional identity for custom functions recorded in the file. Change it when either function's behavior changes. |

The Gaussian is `exp(-(r - k)^2 / (2 * sigma^2))`, where `sigma = params[1]`.
Its default width is `0.7 * length(x) / length(k)`, and its default
integration limit truncates it at radius `k + sqrt(20) * sigma`.

Half-plane bases default to `l=0:4` and use a unit-spaced `x` grid
starting at zero. This grid also defines the radial samples in the
result. To precompute a half-plane basis, use a distinct path.

```julia
basis_path = joinpath(mkpath("data"), "gaussian_half_plane_r512_k128_l4.h5")
basis_data = BASISData(basis_path; geometry=:half_plane, l=0:4, make_images=false)
```

Supply `geometry=:half_plane` to `pbasex` when using this file. A geometry
mismatch raises an error. Existing MATLAB and Julia files without a
geometry field use the quadrant layout. [Forward–backward fitting](usage.md#fit-forwardbackward-asymmetry)

## Use a compact cubic B-spline

Choose `radial_basis=:cubic_bspline` when creating a basis with `pbasex` or
`BASISData`. To select a custom profile instead, omit `radial_basis` and
supply `rBF` and `zIP`.

```julia
using PBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
basis_path = joinpath(mkpath("data"), "cubic_bspline_r512_k128_l4.h5")
result = pbasex(image; basis=basis_path, radial_basis=:cubic_bspline, make_images=false)
```

The basis at radial center `k` is a centered cardinal cubic B-spline
with knot spacing `h = params[1]`. The default is
`h = length(x) / length(k)` pixels, which is 4 px for the default grid and
centers. With $t=(r-k)/h$, the profile is

$$
B_3(t)=\begin{cases}
(4-6|t|^2+3|t|^3)/6 & |t|<1, \\
(2-|t|)^3/6 & 1\le |t|<2, \\
0 & |t|\ge 2.
\end{cases}
$$

The profile and its first two derivatives are continuous at the knots.
Each coefficient affects the radial profile only on
`[k - 2h, k + 2h]`, restricted to physical radii `r >= 0`.

The finite support gives the forward projection exact integration bounds.
Julia maps the support boundaries and interior knots to line-of-sight
coordinates, so the quadrature samples each polynomial piece. Pixels
outside the outer support radius have zero contribution from that basis
function. Pixels inside the projected shell can still receive its
contribution. See [line-of-sight integration](quadrature.md) for the
geometry and the missed-shell example.

Compact support makes radial coefficients more local and lets you control
overlap through `h`. Choose `h` with the center spacing in mind. Narrower
support can leave gaps, while wider support gives more overlap and can
smooth narrow spectral features. Away from the origin, a spline with
`h = 4` px has RMS width `h / sqrt(3)`, about 2.31 px, while the default
Gaussian width on the 512-radius grid is 2.8 px. A comparison using those
defaults measures the combined effect of shape and width. The forward fit
uses a dense SVD, and conditioning depends on the grid and basis overlap.
[Measured spline and Gaussian comparison](matlab-migration.md#circularized-image-with-aligned-polarization)

## Supply custom radial functions

For a custom basis, choose the radial profile and integration limit together.
Passing only `rBF` keeps the Gaussian cutoff. The following small basis
uses a Lorentzian profile and a corresponding finite integration limit.

```julia
using PBASEX

radial(r, k, p) = p[1]^2 / ((r - k)^2 + p[1]^2)
limit(R, k, p) = sqrt(max(0.0, (k + 10 * p[1])^2 - R^2))

basis_path = joinpath(mkpath("data/custom-bases"), "lorentzian-width-2.h5")
custom_basis = BASISData(basis_path; x=0:15, l=[0, 2], params=[2.0],
                         rBF=radial, zIP=limit, cache_key="lorentzian-v1")
```

Reopen a custom basis by path alone to use the saved matrices. To pass
the original functions again, supply the same `cache_key` you used when
creating the file. The HDF5 file stores numeric parameters and metadata.
Keep the Julia definitions of custom functions with your analysis.

## Manage time and memory

A new basis requires many numerical integrals and a dense SVD. Start Julia
with `--threads=auto` to parallelize basis evaluation. Reuse a saved basis
for repeated inversions.

At the default 512 radial samples and 128 Gaussian centers, a quadrant
basis has 262,144 pixels and 384 basis functions. Each dense `Float64`
matrix such as `U` or `Ginv` takes 768 MiB. A half-plane basis with
its default orders has 523,776 pixels and 640 basis functions. Each
such matrix takes about 2.50 GiB.
The SVD and intermediate arrays require additional memory. Increasing the
number of radial centers or angular orders increases the number of columns.

Set `make_images=false` to produce spectra and coefficients without computing
or reading `Ginv`. If a saved basis lacks `Ginv`, a later image-producing run
recomputes the slice matrix from the available radial function. That run
leaves the file as it was.

## HDF5 layout

Let `nx = length(x)`, `nk = length(k)`, `nl = length(l)`, and
`q = length(S)`. The folded detector has `nz = nx` rows for a
quadrant or `nz = 2nx - 1` rows for a half-plane.

| Dataset | Meaning and shape |
| --- | --- |
| `x`, `k`, `l` | Grid, radial centers, and Legendre orders. |
| `radial_basis` | Built-in family identifier. Older files may omit it for Gaussian bases. |
| `geometry` | `quadrant` or `half_plane`. A half-plane file uses `nx * (2nx - 1)` pixels in place of `nx²`, with axial rows ordered from positive to negative. Its slice matrix includes the `1/(4π)` normalization. |
| `params` | Numeric basis parameters. |
| `U`, `S`, `V` | SVD of the projected basis. Sizes are `(nz * nx) × q`, `q`, and `(nk * nl) × q`. |
| `frk` | Radial evaluations with shape `nk × nx`. |
| `Ginv` | Optional central-slice matrix with shape `(nz * nx) × (nk * nl)`. |
| `custom_functions`, `cache_key` | Flags for custom radial and integration functions, and an optional identity string. |

See [MATLAB migration](matlab-migration.md#export-a-julia-basis-for-matlab)
for the CPBASEX export layout.

## Lower-level tools

| Function | Use |
| --- | --- |
| `build_forward_basis(x, k, l, rBF, zIP, params, nothing; geometry)` | Build the forward operator for either geometry and return `U, S, V`. |
| `build_central_slice_basis(x, k, l, rBF, params; geometry)` | Evaluate the stored slice basis. To recover the physical density from `Ginv * vec(result.coefficients)`, divide by `4π` for a quadrant basis. The half-plane matrix includes this factor. `pbasex` handles both conventions. |
| `fold_quadrants(image, column, row; quadrant_filter)` | Sum the selected quadrants for an image or stack. |
| `resize_folded_quadrant(quadrant, nx)` | Crop or zero-pad the folded quadrant to the basis dimensions. |
| `unfold_quadrant(quadrant; center_offset=true, quadrant_count=4)` | Mirror a folded image and divide by the number of quadrants summed. |
| `leg(l, x)` | Evaluate Legendre polynomial `l` at a scalar or array of arguments. |
