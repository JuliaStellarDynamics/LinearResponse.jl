"""
for the radially-biased Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerROI.jl"
include(inputfile)

import CallAResponse
using HDF5

# call the function to construct W matrices
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,d4ψ,FHT,basis,Parameters)

# run the G function calculation
CallAResponse.RunGfunc(ndFdJ,FHT,Parameters)

# construct a grid of frequencies to probe
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

# compute the matrix response at each location in tabomega
tabdet = CallAResponse.RunM(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.2im
bestomg = CallAResponse.FindPole(startingomg,FHT,Parameters,1.e-12)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,FHT,basis,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
