"""
    pBASEX

Julia implementation of the PBASEX algorithm for photoelectron image analysis.

The PBASEX (polar Basis Set Expansion) algorithm analyzes photoelectron angular
distributions by expanding the image in terms of Gaussian-like radial functions
and Legendre polynomials. This allows for the direct extraction of the radial
distribution and angular anisotropy parameters.

This package includes:
- Forward and inverse basis function generators (`findG`, `findGinv`)
- Image processing utilities for folding/unfolding symmetrical images (`foldQuadrant`, `unfoldQuadrant`, `resizeFolded`)
- Memoized Legendre polynomial calculations (`MemoizedLegendre.leg`) (removed)
- Tools for generating basis filenames (`generate_hashed_filename`)
- The main pBASEX algorithm implementation (`pbasex`)

# Key Types and Functions
- `pbasex`: Process photoelectron images using the pBASEX algorithm
- `findG`: Calculate the forward transformation matrix G
- `findGinv`: Calculate the inverse transformation matrix Ginv
- `foldQuadrant`, `unfoldQuadrant`, `resizeFolded`: Image manipulation utilities
- `generate_hashed_filename`: Create unique filenames for basis data
- `PBASEXResult`: Structure holding results from `pbasex`
- `BASISData`: Structure holding basis set definition and precomputed matrices
- `MemoizedLegendre.leg`: Calculate Legendre polynomials (memoized) (removed)
- `clear_cache`, `cache_info`: Manage Legendre polynomial cache (removed)

# Examples
```julia
using pBASEX

# Define basis parameters (example)
basis_params = Dict(
    :x => 0:511,
    :k => (0:511)[2:4:end] .+ 0.5,
    :l => [0, 2, 4]
)

# Create basis data object (will load from file if exists, otherwise calculate)
# Assumes a 'data' directory exists in the package root
basis_data = BASISData(joinpath(pkgdir(pBASEX), "data"); basis_params...)

# Assuming 'my_image' is a loaded 2D image array
# Process the image
result = pbasex(my_image; basis_file=basis_data, options=Dict(:reg => 0.1))

# Extract results
energy = result.E
spectrum = result.IE
beta_params = result.betas
```

# References
- Garcia GA, Nahon L, Powis I. Rev Sci Instrum. 2004;75(11):4989-4996.
  "Two-dimensional charged particle image inversion using a polar basis function expansion"
"""
module pBASEX

# Include sub-modules
include("pbasex_types.jl")
include("math_utils.jl")
#include("memoized_legendre.jl")
include("image_processing.jl")
include("basis_functions.jl")
include("inverse_basis_functions.jl")
include("pbasex.jl")

# Bring necessary modules into scope
using .PBASEXTypes
using .MathUtils 
#using .MemoizedLegendre
using .ImageProcessing
using .BasisFunctions
using .InverseBasisFunctions
using .PBASEX

# --- Public API Exports ---

# Core types
export PBASEXResult, BASISData

# Main algorithm function
export pbasex

# Basis function calculation
export findG, findGinv

# Basis data utilities
export generate_hashed_filename

# Image processing utilities
export foldQuadrant, unfoldQuadrant, resizeFolded

# Legendre polynomial utilities (qualified)
export leg # from MathUtils
#export clear_cache, cache_info

# Convenience function
export greet

# --- Implementation Details ---

# Re-export the specific 'leg' from MathUtils for easier user access
const leg = MathUtils.leg

"""
    greet()

Prints a welcome message for the pBASEX package.
"""
greet() = println("Welcome to pBASEX.jl - A Julia implementation of the PBASEX algorithm")

end # module CPBASEX
