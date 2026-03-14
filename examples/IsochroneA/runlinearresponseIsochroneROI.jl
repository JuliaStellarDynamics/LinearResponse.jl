
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneROIFiducial.jl"
include(inputfile)

import LinearResponse
using HDF5

# compute the Fourier-transformed basis elements
LinearResponse.RunLinearResponse(model,distributionfunction,FHT,basis,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.03im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
