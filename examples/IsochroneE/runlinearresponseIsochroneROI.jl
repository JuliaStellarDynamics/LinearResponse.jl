"""
for the Isochrone radially-unstable model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamIsochroneROIFiducial.jl"
include(inputfile)

import LinearResponse
using HDF5

# combine all steps into one wrapper:
LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.02im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
