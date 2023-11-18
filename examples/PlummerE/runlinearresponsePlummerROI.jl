"""
for the radially-biased Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerROIFiducial.jl"

include(inputfile)

import LinearResponse
using HDF5

# do all linear response steps:
# 1. call the function to construct W matrices
# 2. run the G function calculation
# 3. compute the matrix response coefficients
LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# find a pole by using gradient descent
startingomg = 0.1 + 0.01im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
