# Compatibility entry point for scripts that load the module with include.
if !isdefined(@__MODULE__, :pBASEX)
    include("pBASEX.jl")
end
