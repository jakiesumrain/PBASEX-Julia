# Migrate from MATLAB CPBASEX

Start with the same detector image, center, basis, energy calibration, and
regularization used in MATLAB. A shared basis isolates differences between
the solvers and preserves the original radial sampling.

## Reuse a MATLAB basis in Julia

Pass a CPBASEX HDF5 basis containing `Up` to `pbasex`. Julia reads its grid
and handles the CPBASEX matrix orientations.

```julia
using pBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
matlab_basis = "data/pbasex-matlab/examples/G_r512_k128_l10.h5"
result = pbasex(image; basis=matlab_basis, center=(512, 512),
                alpha=1.0, reg=0.0, make_images=false)
writedlm("energy-spectrum.txt", hcat(result.E, result.IE), '\t')
```

The basis file in this example is a local input. Place your own CPBASEX file
at that path or change `matlab_basis` to match its location. Use the `alpha`
from your MATLAB analysis explicitly, since the two implementations can
have different defaults.

| MATLAB CPBASEX | Julia |
| --- | --- |
| `resizeFolded(foldQuadrant(image,x0,y0,filter),nx)` | Pass the original image with `center=(x0,y0)`, `quadrant_filter=filter`, and the basis file. The same one-based column and row coordinates apply. |
| `loadG(filename, makeImages)` | Pass `basis=filename` and `make_images=makeImages` with a Julia Boolean. |
| `regularization` | `reg` |
| `out.E`, `out.IE`, `out.betas` | `result.E`, `result.IE`, `result.betas` |
| `out.c` | `vec(result.coefficients)` for one image. Julia also exposes the radial-by-angular coefficient matrix. |
| `out.recon`, `out.inv` | `result.reconstructed_image`, `result.inverted_image` |

For an image stack, loop over its slices in Julia. The Julia fit rescales a
subset of quadrants to the four-quadrant intensity scale. Use all four
quadrants when reproducing the validated MATLAB workflow. See the
[usage guide](usage.md#set-the-image-center-and-quadrants) for geometry and
[weighted fits](usage.md#weighted-fits) for how Julia handles weights.

## Export a Julia basis for MATLAB

The following example generates a basis with the same grid, centers,
Gaussian width, and angular orders as the 512-radius comparison. It writes a
CPBASEX-compatible HDF5 file in a separate output directory.

```julia
using pBASEX, HDF5

x = collect(0:511)
k = collect(1.5:4:509.5)
l = collect(0:2:10)
basis = BASISData{Float64}(; x, k, l, params=[2.8])

output = joinpath(mkpath("data/julia-basis"), "G_r512_k128_l10.h5")
h5open(output, "w") do file
    write(file, "x", reshape(Float64.(x), :, 1))
    write(file, "k", reshape(k, 1, :))
    write(file, "l", reshape(Float64.(l), 1, :))
    write(file, "nk", fill(Float64(length(k)), 1, 1))
    write(file, "nl", fill(Float64(length(l)), 1, 1))
    write(file, "Up", basis.U)
    write(file, "S", reshape(basis.S, 1, :))
    write(file, "V", permutedims(basis.V))
    write(file, "frk", basis.frk)
    write(file, "Ginv", permutedims(basis.Ginv))
end
```

With CPBASEX on the MATLAB path, load the file with
`loadG('data/julia-basis/G_r512_k128_l10.h5', 1)`. This export assumes an
unweighted basis and includes the central-slice matrix.

## Validation results

These comparisons used Julia 1.11.3 and MATLAB R2024b with all four quadrants
selected and unweighted fits. Relative L2 error is the norm of the elementwise
difference divided by the norm of the reference output.

| Image and comparison | Result |
| --- | --- |
| Bundled 1023 × 1023 image, both solvers using the same basis | At most `8.0e-13` relative L2 error across the two shared-basis runs for all compared outputs, including spectra, coefficients, anisotropies, and reconstructed images. |
| Bundled image, independently generated bases | Relative L2 error of `1.50e-6` for the energy spectrum and `1.61e-6` for the inverted slice. |
| Separate circularized 1024 × 1024 image, shared MATLAB basis | Relative L2 error of `1.84e-14` for the energy spectrum. |
| Circularized image, independently generated bases | Relative L2 error of `1.18e-6` for the energy spectrum. Both runs placed the main peak at 189 px. |

MATLAB CPBASEX uses trapezoidal integration for the forward basis, while
Julia uses adaptive Gauss–Kronrod integration. Independently generated
bases can therefore differ even when the input parameters agree. The
same-basis comparisons measure solver agreement. Julia's weighted-cache
behavior is covered by package tests. [Validation record in PR #2](https://github.com/jakiesumrain/PBASEX-Julia/pull/2)
