"""
    pBASEX

Polar basis set expansion for inverting velocity-map imaging projections.
Use `pbasex` to fit an image and extract its energy spectrum, angular
anisotropy parameters, and reconstructed images.
"""
module pBASEX

include("math_utils.jl")
include("image_processing.jl")
include("basis_functions.jl")
include("inverse_basis_functions.jl")
include("pbasex_types.jl")
include("inversion.jl")

using .MathUtils
using .ImageProcessing: foldQuadrant, unfoldQuadrant, resizeFolded
using .BasisFunctions: findG
using .InverseBasisFunctions: findGinv
using .PBASEXTypes: PBASEXResult, BASISData, generate_hashed_filename
using .PBASEX: pbasex

export pbasex, PBASEXResult, BASISData
export findG, findGinv, generate_hashed_filename
export foldQuadrant, unfoldQuadrant, resizeFolded, leg, greet

"""
    greet()

Print a welcome message for pBASEX.
"""
greet() = println("Welcome to pBASEX.jl - A Julia implementation of the PBASEX algorithm")

end # module pBASEX
