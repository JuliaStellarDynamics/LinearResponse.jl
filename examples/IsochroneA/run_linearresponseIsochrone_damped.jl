


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

#tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabdet = CallAResponse.RunMIsochrone(inputfile,tabomega,VERBOSE=0)

# for the minimum, go back and compute the mode shape
bestomg = 0.0143 - 0.00141im # for n1max = 10
bestomg = 0.01 - 0.000733im # for n1max = 12
EV,EF,EM = CallAResponse.ComputeModeTablesIsochrone(inputfile,bestomg)
ModeR,ModeShape = CallAResponse.GetModeShape(inputfile,0.01,25.,100,EM)
