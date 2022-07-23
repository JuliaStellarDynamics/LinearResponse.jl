module CallAResponse

# helper functions for constructing basis list, resonance list, prefactors (3d), and writing utilities
include("Utils/Basis.jl")
include("Utils/Resonances.jl")
include("Utils/CMatrix.jl")
include("Utils/IO.jl")

# code to compute the Fourier transform basis functions
include("WMat.jl")

# code to compute the G = n.dFdJ.W.W function
include("GFunc.jl")

# code to compute the Hilbert transformation, including coefficient calculation
include("Xi.jl")

# include code to compute isochrone-specific quantities:
# the isochrone mode is a specific case that we have pre-defined
include("Isochrone/WMatIsochrone.jl")
include("Isochrone/GFuncIsochrone.jl")

end # module
