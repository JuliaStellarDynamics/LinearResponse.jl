
# bring in all parameters from the input file
inputfile = "ModelParamIsochroneROI.jl"
include(inputfile)

import CallAResponse
using HDF5

# call the function to construct W matrices
#CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,FHT,basis,Parameters)

CallAResponse.RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,FHT,basis,Parameters)


# construct a grid of frequencies to probe
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)


# compute the matrix response at each location
#tabdet = CallAResponse.RunM(tabomega,dψ,d2ψ,FHT,basis,Parameters)



bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,dψ,d2ψ,FHT,basis,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
