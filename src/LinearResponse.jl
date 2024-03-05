module LinearResponse

# bring in the external dependencies
using AstroBasis
using FiniteHilbertTransform
using HDF5
using LinearAlgebra # Access to Symmetric
using OrbitalElements

# Default values
const defaultΩ₀ = 1.
const defaultrmin = 1.0e-6
const defaultrmax = 1.0e4

# structure to hold all parameters
include("Utils/ParameterStructure.jl")


# helper functions for constructing basis list, resonance list, prefactors (3d), writing utilities, and frequency tables
include("Utils/Resonances.jl")
include("Utils/CMatrix.jl")
include("Utils/IO.jl")
include("Utils/Compatibilities.jl")
include("Utils/OmegaGrid.jl")

# code to compute the Fourier transform basis functions
include("WMat.jl")

# code to compute the G = n.dFdJ.W.W function
include("GFunc.jl")

# code to compute the Hilbert transformation, including coefficient calculation
include("Xi.jl")

# Code to compute the response matrix
include("ResponseMatrix.jl")

# Codes to extract different informations
include("ModeComputation/Determinant.jl")
include("ModeComputation/Matrices.jl")
include("ModeComputation/Mode.jl")
include("ModeComputation/FindPole.jl")

# include code to compute isochrone-specific quantities:
# the isochrone mode is a specific case that we have pre-defined
include("Isochrone/WMatIsochrone.jl")

include("TimeResponse.jl")

end # module
