# PBASEX-Julia

`PBASEX` uses polar basis set expansion to invert velocity-map imaging (VMI)
projections. It fits a detector image and returns the energy spectrum,
angular anisotropy coefficients, and optionally the reconstructed projection
and Abel-inverted central slice.

PBASEX supports two detector geometries. Quadrant mode uses both
reflection symmetries and fits even Legendre orders. Half-plane mode
preserves forward–backward asymmetry along an in-plane physical axis and
fits both even and odd orders.

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

Pass a basis filename to every inversion. PBASEX creates the file at a new
path and reuses it on later runs. To change the grid or radial basis,
choose a new filename. You can also pass the path to a Julia or
MATLAB basis you already have. [Choose or reuse a basis](docs/bases.md#select-a-basis)

The default energy axis is radius squared in pixel². Pass `alpha` to convert
it to your calibrated energy units. For input in detected-electron counts,
both geometries return a spectrum in electrons per unit energy, a fitted
projection in electrons per detector pixel, and a central slice of the 3D
density in electrons per cubic pixel. [Output units and peak yields](docs/usage.md#physical-units-and-peak-yields)

For your own images, align the polarization axis with the vertical image
axis before fitting. The Gaussian basis is the default. Select
`radial_basis=:cubic_bspline` when building a compact spline basis.
[Image alignment](docs/usage.md#align-the-polarization-axis) · [Spline basis](docs/bases.md#use-a-compact-cubic-b-spline)

For a cylindrically symmetric image with distinct upper and lower
hemispheres, orient the positive physical axis toward the top and use
half-plane mode with its own basis file.

```julia
half_plane_basis = joinpath(mkpath("data"), "gaussian_half_plane_r512_k128_l4.h5")
result = pbasex(image; basis=half_plane_basis, geometry=:half_plane, l=0:4,
                make_images=false)
beta_orders = result.legendre_orders[2:end]
```

[Half-plane fitting and sign convention](docs/usage.md#fit-forwardbackward-asymmetry)

Run the bundled image reconstruction example from the repository root.
The first script prepares `basis.h5`, and the second uses it to invert the
bundled detector image and write `out.h5`.

```sh
julia --project=. --threads=auto examples/generate_basis.jl
julia --project=. --threads=auto examples/invert_detector_image.jl
```

## Documentation

- [Run an inversion](docs/usage.md) covers image coordinates, quadrant selection,
  regularization, and output fields.
- [Build and save a basis](docs/bases.md) covers radial and angular settings,
  file reuse, memory requirements, and the HDF5 layout.
- [Line-of-sight integration](docs/quadrature.md) explains how missed
  quadrature samples can produce spurious low-energy peaks.
- [Migrate from MATLAB CPBASEX](docs/matlab-migration.md) covers using an existing
  MATLAB basis, exporting a compatible basis from Julia, and the measured parity.
- [Original algorithm derivation](docs/pbasex.pdf) describes the Abel transform
  and even-order polar basis expansion.
- [From detected electrons to pBASEX error bars](docs/uncertainty-derivation.md)
  derives the spectrum, window yields and beta uncertainties from the original
  pBASEX equations, with worked examples of pixel counting noise and shot resampling.
- [Forward–backward asymmetry](docs/forward-backward-asymmetry.md) derives
  the implemented half-plane geometry and its odd Legendre terms.
- [Counting-error estimation](docs/counting-errors.md) uses `pbasex_counts` to estimate uncertainties
  for the spectrum and anisotropies at each energy, and for yields and
  anisotropies integrated over a peak. The guide also shows how to report
  [a spectrum and betas in uniform energy bins](docs/counting-errors.md#yield-preserving-spectrum-and-betas-on-a-uniform-energy-grid).
- [How shot-resolved pBASEX works](docs/shot-uncertainties.md)
  uses `pbasex_shots` to bootstrap arbitrary energy intervals from corrected event checkpoints.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Citation and license

If you use this software in research, cite García, Nahon, and Powis (2004),
*Review of Scientific Instruments* 75, 4989–4996.
[Original pBASEX paper](https://doi.org/10.1063/1.1807578)

PBASEX-Julia builds on Elio Gilbert Champenois's CPBASEX implementation in
MATLAB and Python. The Julia package adds compact cubic B-spline radial bases,
half-plane inversion with odd Legendre orders for forward–backward asymmetry,
and uncertainty estimates for electron spectra and angular anisotropies.
It propagates Poisson counting errors from detector histograms and estimates
shot-to-shot uncertainty by bootstrapping corrected event checkpoints.
[CPBASEX repository](https://github.com/e-champenois/CPBASEX)

PBASEX-Julia is distributed under the MIT license. [License](LICENSE)
