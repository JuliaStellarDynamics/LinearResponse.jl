module CallAResponse

# bring in the external dependencies
import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
using HDF5


# helper functions for constructing basis list, resonance list, prefactors (3d), writing utilities, and frequency tables
include("Utils/Basis.jl")
include("Utils/Resonances.jl")
include("Utils/CMatrix.jl")
include("Utils/IO.jl")
include("Utils/OmegaGrid.jl")

# code to compute the Fourier transform basis functions
include("WMat.jl")

# code to compute the G = n.dFdJ.W.W function
include("GFunc.jl")

# code to compute the Hilbert transformation, including coefficient calculation
include("Xi.jl")

# code to compute the shape of the mode
include("Mode.jl")

# include code to compute isochrone-specific quantities:
# the isochrone mode is a specific case that we have pre-defined
include("Isochrone/WMatIsochrone.jl")
include("Isochrone/GFuncIsochrone.jl")
include("Isochrone/XiIsochrone.jl")

end # module
