


import CallAResponse
using HDF5

# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochrone_roi.jl"

# compute the Fourier-transformed basis elements
CallAResponse.RunWmatIsochrone(inputfile)

# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(inputfile)

include(inputfile)
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunMIsochrone(inputfile,tabomega,VERBOSE=1)

#=
# compute the determinants with a gradient descent to find the exact frequency for ROI
bestomg = CallAResponse.FindZeroCrossing(inputfile,0.00,0.03,NITER=16,VERBOSE=1)

# for n1max=3
#bestomg = 0.0 + 0.02271406012170436im

# for n1max=10
#bestomg = 0.0 + 0.022914332993273286im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(inputfile,bestomg)
ModeR,ModeShape = CallAResponse.GetModeShape(inputfile,0.01,15.,100,EM)

=#
