


import CallAResponse
using HDF5

# this input file is set up to mimic the Fouvry & Prunet (2022) damped mode result exactly
inputfile = "ModelParamIsochrone_damped.jl"

# compute the Fourier-transformed basis elements
CallAResponse.RunWmatIsochrone(inputfile)

# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(inputfile)

# need this to get the parameters...
include(inputfile)

tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunMIsochrone(inputfile,tabomega,VERBOSE=2)
