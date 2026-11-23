# PBASEX-Julia

> ⚠️ **Note**: This README was generated with AI assistance to provide comprehensive documentation based on the source code.

A Julia implementation of the PBASEX (Polar Basis Set Expansion) algorithm for analyzing 2D photoelectron/photoion images from velocity map imaging (VMI) experiments.

## Overview

PBASEX decomposes 2D projection images using:
- **Radial basis functions**: Gaussians centered at different radii
- **Angular basis functions**: Legendre polynomials P_l(cos θ)

This allows reconstruction of the original 3D distribution from 2D projections, extracting:
- Energy spectrum I(E) - radial intensity distribution
- Angular anisotropy parameters β_l(E) - describes angular distribution

Based on: Garcia, G. A., Nahon, L., & Powis, I. (2004). Two-dimensional charged particle image inversion using a polar basis function expansion. *Review of Scientific Instruments*, *75*(11), 4989-4996. [doi:10.1063/1.1807578](https://doi.org/10.1063/1.1807578)

**This is a refactoring and Julia port of the original MATLAB/Python implementation**: [github.com/e-champenois/CPBASEX](https://github.com/e-champenois/CPBASEX)

## Features

- **Automatic basis caching**: Computationally expensive basis matrices are cached as HDF5 files
- **Multi-threading support**: Parallel computation of basis integrals
- **Flexible basis functions**: Customizable radial and angular basis functions
- **Image preprocessing**: Automatic quadrant folding and resizing

## Installation

Clone the repository and instantiate its project with Julia 1.10 or later in
the Julia 1 series.

```bash
git clone https://github.com/jakiesumrain/PBASEX-Julia.git
cd PBASEX-Julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Dependencies

The project declares its dependencies in `Project.toml`.
- **HDF5.jl**: Basis data caching
- **QuadGK.jl**: Numerical integration
- **ProgressMeter.jl**: Progress bars during computation
- **LinearAlgebra** (stdlib)
- **Serialization, SHA** (stdlib): Hash-based caching
- **DelimitedFiles**: Text image input in the examples

## Quick Start

### Basic Example

```julia
# Start Julia with julia --project=. from the repository root
using pBASEX
using DelimitedFiles, HDF5

# Define the basis function parameters
l_values = collect(0:2:4)  # Legendre polynomial orders: 0, 2, 4

# Load a square floating-point image
img = readdlm("examples/image.txt", ' ', Float64)
println(size(img))

# Run PBASEX analysis
result = pbasex(img, l=l_values, basis="basis.h5")

# Access results
energy = result.E                       # Energy scale (proportional to r²)
intensity = result.IE                   # I(E) spectrum
beta_2 = result.betas[:, 1]            # β₂(E) anisotropy parameter
beta_4 = result.betas[:, 2]            # β₄(E) anisotropy parameter
recon_img = result.reconstructed_image # Reconstructed projection
inv_img = result.inverted_image        # Abel inverted image (central slice)

# Save results to HDF5 file
h5open("output.h5", "w") do file
    write(file, "E", result.E)
    write(file, "IE", result.IE)
    write(file, "betas", result.betas)
    write(file, "inv_img", result.inverted_image)
end

# Note: HDF5 files saved in Python may have transposed matrices when read in Julia (and vice versa)
```

### Running Examples

**With multi-threading** (recommended):
```bash
julia --project=. --threads=auto examples/sample_pbasex.jl
```

**Single-threaded**:
```bash
julia --project=. examples/sample_pbasex.jl
```

Run the test suite from the repository root.

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Usage Guide

### The `pbasex` Function

**Signature:**
```julia
pbasex(image::AbstractMatrix{T};
       basis::Union{AbstractString, Nothing} = nothing,
       options...) where T<:AbstractFloat
```

**Required Input:**
- `image`: Square 2D floating-point array. Use `Float64` for working with
  full-precision bases or `Float32` to reduce matrix memory.

**Returns:** `PBASEXResult` struct with fields:
- `E::Vector`: Energy scale (α × x²)
- `IE::Vector`: Intensity spectrum I(E)
- `betas::Matrix`: Anisotropy parameters β_l(E), columns correspond to l>0 orders
- `reconstructed_image::Matrix`: Reconstructed 2D projection (if `make_images=true`)
- `inverted_image::Matrix`: Abel inverted central slice (if `make_images=true`)
- `residuals::Matrix`: Original minus reconstructed image, when their size and center match
- `coefficients::Matrix`: Raw basis expansion coefficients [nk × nl]

### Key Parameters

#### Basis Control

**`basis`** - Controls basis caching:
- `nothing` (default): Uses the repository's `data/` directory.
- `"directory/"`: Uses an existing directory. Its cache key includes the
  grid, basis parameters, and processed pixel weights. For custom radial or
  integration functions, supply a `cache_key` that identifies their behavior.
  Without one, the custom basis is computed in memory.
- `"file.h5"`: Creates the named file or loads its existing basis and grid.
  Explicitly supplied parameters must match the saved definition.
  MATLAB CPBASEX files with an `Up` dataset can be used directly.

#### Basis Definition

**`x`** - Radial pixel grid (default: `0:511`)
```julia
x = collect(0:511)  # Standard 512-pixel radius
```

**`k`** - Radial basis function centers (default: `x[2:4:end] .+ 0.5`)
```julia
k = collect(0:2:511) .+ 0.5  # Denser radial sampling
```

**`l`** - Even Legendre polynomial orders starting with zero (default: `[0, 2, 4]`)
```julia
l = [0, 2, 4, 6]  # Add β₆ parameter
```

**`params`** - Gaussian width parameter (default: `[0.7 * length(x) / length(k)]`)
```julia
params = [1.5]  # Wider Gaussian basis functions
```

**`rBF`** - Radial basis function (default: Gaussian)
```julia
# Default: Gaussian
rBF = (r, k, p) -> exp(-(r - k)^2 / (2 * p[1]^2))

# Custom: Exponential
rBF = (r, k, p) -> exp(-abs(r - k) / p[1])
```

**`zIP`** - Integration upper limit (default: based on Gaussian properties)
```julia
# Default for Gaussian
zIP = (R, k, p) -> sqrt(max(0.0, (sqrt(20.0)*p[1] + k)^2 - R^2))
```

#### Processing Control

**`reg`** - Tikhonov regularization parameter (default: `0.0`)
```julia
reg = 0.01  # Stabilize noisy images (typical range: 0.001-0.1)
```

**`make_images`** - Reconstruct images (default: `true`)
```julia
make_images = false  # Skip image reconstruction for speed
```

This also skips building or loading the dense central-slice matrix `Ginv`.

**`center`** - Image center in (column, row) (default: geometric center)
```julia
center = (257, 257)  # Manually specify center
```

**`quadrant_filter`** - Quadrant selection (default: `[1, 1, 1, 1]`)
```julia
quadrant_filter = [1, 1, 1, 1]  # Use all four quadrants
quadrant_filter = [1, 0, 0, 0]  # Use only top-right quadrant
# MATLAB CPBASEX order: [Top-Right, Top-Left, Bottom-Left, Bottom-Right]
```

Select at least one quadrant. The inversion normalizes the selected quadrants
to the intensity scale of a four-quadrant sum.

**`alpha`** - Energy scaling factor (default: `1.0`)
```julia
alpha = 2.5  # Scale energy axis: E = alpha * x²
```

**`weights`** - Pixel weighting matrix (default: `nothing`)
```julia
# Must match original image size before folding/resizing
weights = ones(512, 512)  # Uniform weighting
```
Pass the same weight image when reusing a named weighted Julia basis. Files
written by older Julia versions without weight metadata are treated as
unweighted. Regenerate a legacy weighted file before reusing it.

### Complete Example

```julia
using pBASEX
using DelimitedFiles

# Load experimental image
img = readdlm("my_image.txt", Float64)

# Configure analysis
result = pbasex(img,
    # Basis caching
    basis = "cache/my_basis.h5",

    # Basis definition
    x = collect(0:511),
    k = collect(0:2:511) .+ 0.5,  # Denser sampling
    l = [0, 2, 4, 6],             # Include β₆

    # Processing
    reg = 0.05,                    # Regularization for noisy data
    center = (256, 256),           # Off-center by 1 pixel
    alpha = 2.0,                   # Energy scaling
    make_images = true             # Generate reconstructed images
)

# Save results
using HDF5
h5open("results.h5", "w") do file
    write(file, "energy", result.E)
    write(file, "intensity", result.IE)
    write(file, "betas", result.betas)
    write(file, "inverted_image", result.inverted_image)
end

# Plot results
using Plots
plot(result.E, result.IE, label="I(E)", xlabel="Energy", ylabel="Intensity")
plot!(result.E, result.betas[:, 1], label="β₂(E)")
```

## Advanced Usage

### Pre-computing Basis Sets

For repeated analyses with the same basis configuration, pre-compute and cache the basis:

```julia
using pBASEX

# Define the basis function parameters
l_values = collect(0:2:4)  # Legendre polynomial orders
x_values = collect(0:511)  # Radial grid points

# Compute and save basis to a specific file
basis_data = BASISData{Float64}("./basis.h5", x=x_values, l=l_values)

# Or save to a directory (hash-based filename)
# basis_data = BASISData{Float64}("./", x=x_values, l=l_values)

# Later, this basis will be loaded from cache automatically
result = pbasex(img, basis="./basis.h5")
```

### Custom Radial Basis Functions

```julia
# Example: Lorentzian basis instead of Gaussian
my_rBF = (r, k, p) -> p[1]^2 / ((r - k)^2 + p[1]^2)
my_zIP = (R, k, p) -> sqrt(max(0.0, (10.0*p[1] + k)^2 - R^2))

result = pbasex(img,
                basis = "lorentzian_basis.h5",
                cache_key = "lorentzian-v1",
                rBF = my_rBF,
                zIP = my_zIP,
                params = [2.0])
```

### Weighted Inversion

Apply pixel-specific weights (e.g., for noise statistics):

```julia
# Create weight matrix (same size as original image)
noise_std = estimate_noise(img)
weights = 1.0 ./ noise_std

result = pbasex(img,
                basis = "weighted_basis.h5",
                weights = weights,
                reg = 0.01)
```

## Performance Optimization

### Multi-threading

Enable Julia threads for parallel basis computation:
```bash
# Automatic thread detection
julia --project=. --threads=auto script.jl

# Specify thread count
julia --project=. --threads=8 script.jl
```

### Basis Caching Strategy

First run with new basis parameters is slow (minutes to hours depending on grid size):
```julia
# First run: integrates and factorizes the basis
result1 = pbasex(img1, basis="shared_basis.h5", x=collect(0:511), l=[0,2,4])

# Subsequent runs: loads from cache (~seconds)
result2 = pbasex(img2, basis="shared_basis.h5")
result3 = pbasex(img3, basis="shared_basis.h5")
```

### Memory Considerations

For large images or dense basis sets:
- Use `make_images=false` if only spectra and anisotropies are needed.
- Each dense matrix such as `U` or `Ginv` uses roughly
  `sizeof(T) * length(x)^2 * length(k) * length(l)` bytes for a tall basis.

## Understanding Results

### Energy Scale
```julia
E = result.E  # E[i] = α × x[i]²
```
Units depend on your image calibration and the `alpha` parameter.

### Intensity Spectrum
```julia
IE = result.IE  # I(E) - 1D energy distribution
```
This is the l=0 (isotropic) component.

### Beta Parameters
```julia
betas = result.betas  # [n_energies × n_betas]
# betas[:, 1] = β₂(E)  - if l includes 2
# betas[:, 2] = β₄(E)  - if l includes 4
# etc.
```

Angular distribution: I(E, θ) = I(E) × [1 + Σ β_l(E) × P_l(cos θ)]

### Images
```julia
reconstructed_image  # 2D projection from fit coefficients
inverted_image       # Abel inverted central slice
residuals            # Original - reconstructed (fit quality)
```

## Troubleshooting

### "Input image must be square"
Ensure your image dimensions are equal: `size(img) = (n, n)`

### "Parent directory for specified basis_file does not exist"
Create the directory first: `mkpath("path/to/cache")`

### Poor fit quality (large residuals)
- Increase regularization: `reg = 0.01` to `0.1`
- Adjust basis density: more k points for finer features
- Verify image center: `center=(col, row)`

### Slow computation
- Enable threading: `julia --threads=auto`
- Reduce basis size: use sparser k grid
- Use cached basis files

### Type errors
Convert integer images before inversion, for example with `img = Float64.(img)`.

## HDF5 File Format

### Basis Cache Files
```
BASISData_<hash>.h5:
├── x          # Radial grid [Int]
├── k          # Basis centers
├── l          # Legendre orders [Int]
├── params     # Basis parameters
├── weights    # Processed pixel weights, when supplied
├── U          # Left singular vectors
├── S          # Singular values
├── V          # Right singular vectors
├── frk        # Cached radial basis evaluations
├── Ginv       # Central-slice matrix, when saved
├── custom_functions
└── cache_key  # Custom basis identity, when supplied
```

### MATLAB CPBASEX bases

Pass a CPBASEX HDF5 basis directly to `pbasex`. For an existing file, Julia
reads the coordinate grid and spectral factors and handles the MATLAB
orientations of `Up`, `V`, `frk`, and `Ginv`.

```julia
result = pbasex(img; basis="G_r512_k128_l10.h5")
```

### Results Files (example)
```julia
h5open("output.h5", "w") do file
    write(file, "E", result.E)
    write(file, "IE", result.IE)
    write(file, "betas", result.betas)
    write(file, "inverted_image", result.inverted_image)
end
```

**Note**: HDF5 matrices written by Python may appear transposed in Julia and vice versa.

## Package Structure

```
PBASEX-Julia/
├── Project.toml                   # Package metadata and dependencies
├── src/
│   ├── pBASEX.jl                   # Package entry point
│   ├── inversion.jl                # pbasex() function, preprocessing
│   ├── pbasex_types.jl             # BASISData, PBASEXResult definitions
│   ├── basis_functions.jl          # findG() - forward transform
│   ├── inverse_basis_functions.jl  # findGinv() - inverse transform
│   ├── image_processing.jl         # Folding, resizing utilities
│   └── math_utils.jl               # Legendre polynomials, meshgrid
├── examples/
│   ├── sample_pbasex.jl           # Complete analysis example
│   ├── sample_basis.jl            # Basis pre-computation example
│   └── image.txt                  # Sample data
├── docs/                          # A brief introduction to pbasex method
├── test/
│   └── runtests.jl                 # Package tests
├── LICENSE
└── README.md
```

## Citation

If you use this code, please cite the original PBASEX paper:

```bibtex
@article{garcia2004pbasex,
  title={Two-dimensional charged particle image inversion using a polar basis function expansion},
  author={Garcia, Gustavo A and Nahon, Laurent and Powis, Ivan},
  journal={Review of Scientific Instruments},
  volume={75},
  number={11},
  pages={4989--4996},
  year={2004},
  doi={10.1063/1.1807578}
}
```

## License

MIT License - See [LICENSE](LICENSE) file for details.

This project is a Julia refactoring of the original CPBASEX code:
- Original MATLAB/Python implementation: Copyright (c) 2017 Elio Gilbert Champenois
- Julia refactoring: Copyright (c) 2025 TAO Jianfei

## Author

**TAO Jianfei**
Email: jakiesumrain@163.com

## Contributing

Contributions are welcome! Please submit issues or pull requests on GitHub.
