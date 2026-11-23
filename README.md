# PBASEX-Julia

`pBASEX` uses polar basis set expansion to invert velocity-map imaging (VMI)
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
using pBASEX, DelimitedFiles

image = readdlm("examples/image.txt", Float64)
result = pbasex(image; make_images=false)
writedlm("energy-spectrum.txt", hcat(result.E, result.IE), '\t')
```

The first run builds a basis and caches it in `data/`. Later runs can reuse
that cache. The default energy axis is radius squared in pixel². Pass
`alpha` to convert it to your calibrated energy units.

To run the bundled example with image reconstruction and HDF5 output, use
`julia --project=. --threads=auto examples/sample_pbasex.jl` from the checkout.
It uses `basis.h5` and writes `out.h5` in the working directory.

## Documentation

- [Run an inversion](docs/usage.md) covers image coordinates, quadrant selection,
  regularization, weights, and output fields.
- [Build and cache a basis](docs/bases.md) covers radial and angular settings,
  cache identity, memory requirements, and the HDF5 layout.
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
implementation. See the [MIT license](LICENSE) for copyright details and the
[CPBASEX repository](https://github.com/e-champenois/CPBASEX) for the upstream code.
