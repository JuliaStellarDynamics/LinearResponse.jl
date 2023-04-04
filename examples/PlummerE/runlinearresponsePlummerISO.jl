"""
for the isotropic Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerISO.jl"
include(inputfile)

import LinearResponse
using HDF5

# call the function to construct W matrices
LinearResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,d4ψ,FHT,basis,Parameters)

# run the G function calculation
LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# construct a grid of frequencies to probe
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

# compute the matrix response at each location in tabomega
tabdet = LinearResponse.RunM(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.1 - 0.001im
#bestomg = LinearResponse.FindPole(startingomg,FHT,Parameters,1.e-12)

#bestomg = 0.0 + 0.02271406012170436im
#println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
#EV,EF,EM = LinearResponse.ComputeModeTables(bestomg,FHT,basis,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
#ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
