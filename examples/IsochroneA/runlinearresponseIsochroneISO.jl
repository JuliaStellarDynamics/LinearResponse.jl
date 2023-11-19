
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneISOFiducial.jl"
include(inputfile)

import LinearResponse
using HDF5

# compute the Fourier-transformed basis elements
LinearResponse.RunWmatIsochrone(FHT,bc,M,G,basis,Parameters)

# compute the G(u) functions
LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# compute the matrix response at each location in tabomega
tabdet = LinearResponse.RunAXi(FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.014 + 0.0001im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
