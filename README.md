# PBASEX-Julia

`PBASEX` uses polar basis set expansion to invert velocity-map imaging (VMI)
projections. It fits a square detector image and returns the energy spectrum,
angular anisotropy coefficients, and optionally the reconstructed projection
and Abel-inverted central slice.

## Quick start

Clone the repository, then activate it with Julia 1.10 or later in the Julia 1
series.

```sh
git clone https://github.com/jakiesumrain/PBASEX-Julia.git
cd PBASEX-Julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=auto
```

At the Julia prompt, extract the spectrum from the bundled image and write
two columns to a text file.

```julia
using PBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
basis_path = joinpath(mkpath("data"), "gaussian_r512_k128_l4.h5")
result = pbasex(image; basis=basis_path, make_images=false)
writedlm("energy-spectrum.txt", hcat(result.E, result.IE), '\t')
```

Choose a basis filename and pass its path to every inversion. The first run
creates that file, and later runs reuse it. The default energy axis is radius
squared in pixel². Pass
`alpha` to convert it to your calibrated energy units.

For your own images, align the polarization axis with the vertical image
axis before fitting. The Gaussian basis is the default. Select
`radial_basis=:cubic_bspline` when building a compact spline basis.
[Image alignment](docs/usage.md#align-the-polarization-axis) · [Spline basis](docs/bases.md#use-a-compact-cubic-b-spline)

Run the bundled image reconstruction example from the repository root.
The first script prepares `basis.h5`, and the second uses it to invert the
bundled detector image and write `out.h5`.

```sh
julia --project=. --threads=auto examples/generate_basis.jl
julia --project=. --threads=auto examples/invert_detector_image.jl
```

## Documentation

- [Run an inversion](docs/usage.md) covers image coordinates, quadrant selection,
  regularization, weights, and output fields.
- [Build and save a basis](docs/bases.md) covers radial and angular settings,
  file reuse, memory requirements, and the HDF5 layout.
- [Line-of-sight integration](docs/quadrature.md) explains how missed
  quadrature samples can produce spurious low-energy peaks.
- [Migrate from MATLAB CPBASEX](docs/matlab-migration.md) covers using an existing
  MATLAB basis, exporting a compatible basis from Julia, and the measured parity.
- [Algorithm derivation](docs/pbasex.pdf) describes the Abel transform and polar
  basis expansion.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Citation and license

If you use this software in research, cite García, Nahon, and Powis (2004),
*Review of Scientific Instruments* 75, 4989–4996.
[Original pBASEX paper](https://doi.org/10.1063/1.1807578)

PBASEX-Julia is a Julia port of Elio Gilbert Champenois's MATLAB and Python
implementation, distributed under the MIT license. [License](LICENSE) · [Upstream CPBASEX repository](https://github.com/e-champenois/CPBASEX)
