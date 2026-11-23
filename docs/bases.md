# Build and cache a basis

A basis combines radial profiles with Legendre polynomials. The forward
matrix projects each basis function onto the detector. The solver stores
its singular value decomposition and uses the same fit coefficients to
evaluate the radial profiles and central slice.

## Select a basis

| `basis` argument to `pbasex` | Behavior |
| --- | --- |
| `nothing` | Use the package's `data/` directory. |
| Existing directory | Reuse or create a file keyed by the grid, basis parameters, and processed weights. Custom functions also require a stable `cache_key` for reuse. |
| Existing HDF5 file | Read the grid and factors from that file. Explicitly supplied settings must agree with the file's metadata. Julia and MATLAB CPBASEX layouts are supported. |
| New HDF5 path | Compute and save the basis at that path. The parent directory must exist. |

Precompute a basis once and use its path across images. The ordinary
`BASISData(path)` constructor uses `Float64`. Use `BASISData{Float32}` for
single-precision storage.

```julia
using pBASEX

basis_path = joinpath(mkpath("data"), "even-orders-0-4.h5")
basis_data = BASISData(basis_path; x=0:511, l=[0, 2, 4])
# Reuse the path after loading an image.
# result = pbasex(image; basis=basis_path, make_images=false)
```

Basis reuse is independent of `alpha` and `reg`. For named weighted bases,
pass the same original weights and preprocessing options when calling
`pbasex`. If you precompute a weighted basis with `BASISData`, pass the
processed weight quadrant of size `length(x) × length(x)`. Legacy Julia files
without weight metadata load as unweighted. Regenerate any legacy basis
originally built with weights.

## Set the basis definition

| Keyword | Default for a new basis | Meaning |
| --- | --- | --- |
| `x` | `0:511` | Integer radial grid. Use `0:rmax` for a detector sampled at unit pixel spacing. |
| `k` | `x[2:4:end] .+ 0.5` | Radial centers. |
| `l` | `[0, 2, 4]` | Distinct nonnegative even Legendre orders, starting with zero for image inversion. |
| `params` | `[0.7 * length(x) / length(k)]` | Parameters for both basis functions. For the default Gaussian, the first value is the width in pixels. |
| `rBF` | Gaussian | Radial function with signature `rBF(r, k, params)`. |
| `zIP` | Gaussian cutoff | Upper limit of the line-of-sight integral, with signature `zIP(R, k, params)`. Here `R` is the projected radius. |
| `cache_key` | `nothing` | String identifying a custom radial function and integration limit. Change it when either function's behavior changes. |

The Gaussian is `exp(-(r - k)^2 / (2 * sigma^2))`, where `sigma = params[1]`.
The default integration limit truncates it at radius `k + sqrt(20) * sigma`.

For a custom basis, choose the radial profile and integration limit together.
The following small basis uses a Lorentzian profile.

```julia
using pBASEX

radial(r, k, p) = p[1]^2 / ((r - k)^2 + p[1]^2)
limit(R, k, p) = sqrt(max(0.0, (k + 10 * p[1])^2 - R^2))

cache = mkpath("data/custom-bases")
custom_basis = BASISData(cache; x=0:15, l=[0, 2], params=[2.0],
                         rBF=radial, zIP=limit, cache_key="lorentzian-v1")
```

A custom basis without a key is computed in memory when you pass a directory.
An explicit file can still store the computed matrices. Reopen that file
by path alone to use the saved matrices, or supply the original custom
functions together with their matching `cache_key`. The HDF5 file stores
numeric parameters and metadata. Keep the Julia definitions of custom
functions with your analysis.

## Manage time and memory

A new basis requires many numerical integrals and a dense SVD. Start Julia
with `--threads=auto` to parallelize basis evaluation. Reuse a saved basis
for repeated inversions.

A default-size basis has 262,144 quadrant pixels and 384 basis functions.
Each dense `Float64` matrix such as `U` or `Ginv` then takes 768 MiB.
The SVD and intermediate arrays require additional memory. Increasing the
number of radial centers or angular orders increases the number of columns.

Set `make_images=false` to produce spectra and coefficients without computing
or reading `Ginv`. If a saved basis lacks `Ginv`, a later image-producing run
recomputes the slice matrix from the available radial function. That run
leaves the file as it was.

## HDF5 layout

For a basis with `nx = length(x)`, `nk = length(k)`, `nl = length(l)`, and
`q = length(S)`, Julia's native cache stores these datasets.

| Dataset | Meaning and shape |
| --- | --- |
| `x`, `k`, `l` | Grid, radial centers, and Legendre orders. |
| `params` | Numeric basis parameters. |
| `weights` | Optional processed weight matrix of size `nx × nx`. |
| `U`, `S`, `V` | SVD of the projected basis. Sizes are `nx² × q`, `q`, and `(nk * nl) × q`. |
| `frk` | Radial evaluations with shape `nk × nx`. |
| `Ginv` | Optional central-slice matrix with shape `nx² × (nk * nl)`. |
| `custom_functions`, `cache_key` | Flags for custom radial and integration functions, and an optional identity string. |

See [MATLAB migration](matlab-migration.md#export-a-julia-basis-for-matlab)
for the CPBASEX export layout.

## Lower-level tools

| Function | Use |
| --- | --- |
| `findG(x, k, l, rBF, zIP, params, weights)` | Build the forward operator and return its factors `U, S, V`. Pass `nothing` for an unweighted basis. |
| `findGinv(x, k, l, rBF, params)` | Evaluate the central-slice basis. Multiply this matrix by the fitted coefficient vector. |
| `foldQuadrant(image, column, row; quadrant_filter)` | Sum the selected quadrants for an image or stack. |
| `resizeFolded(quadrant, nx)` | Crop or zero-pad the folded quadrant to the basis dimensions. |
| `unfoldQuadrant(quadrant; center_offset=true, quadrant_count=4)` | Mirror a folded image and divide by the number of quadrants summed. |
| `leg(l, x)` | Evaluate Legendre polynomial `l` at a scalar or array of arguments. |
| `generate_hashed_filename(x, k, l)` | Generate a hash from the three vectors. Automatic basis caches use additional metadata in their keys. |
