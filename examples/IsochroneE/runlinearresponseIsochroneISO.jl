"""
for the Isochrone isotropic model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamIsochroneISO.jl"
inputfile = "ModelParamIsochroneISOFiducial.jl"

include(inputfile)

import LinearResponse
using HDF5

# combine all steps into one wrapper:
LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
#tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.02 - 0.002im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg, with determinant $detval")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)
println(EM)

modeRmin = 0.01
modeRmax = 10.0
nmode = 400
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
