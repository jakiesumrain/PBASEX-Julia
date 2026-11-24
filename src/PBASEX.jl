"""
    PBASEX

Polar basis set expansion for inverting velocity-map imaging projections.
Use `pbasex` to fit an image and extract its energy spectrum, angular
anisotropy parameters, and reconstructed images.
"""
module PBASEX

include("math_utils.jl")
include("image_processing.jl")
include("forward_basis.jl")
include("central_slice_basis.jl")
include("pbasex_types.jl")
include("inversion.jl")

using .MathUtils
using .ImageProcessing: fold_quadrants, unfold_quadrant, resize_folded_quadrant
using .BasisFunctions: build_forward_basis
using .InverseBasisFunctions: build_central_slice_basis
using .PBASEXTypes: PBASEXResult, BASISData
using .Inversion: pbasex

export pbasex, PBASEXResult, BASISData
export build_forward_basis, build_central_slice_basis
export fold_quadrants, unfold_quadrant, resize_folded_quadrant, leg, greet

"""
    greet()

Print a welcome message for PBASEX.
"""
greet() = println("Welcome to PBASEX.jl - A Julia implementation of the PBASEX algorithm")

end # module PBASEX
