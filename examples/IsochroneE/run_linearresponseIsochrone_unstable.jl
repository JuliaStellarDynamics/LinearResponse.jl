


import CallAResponse
using HDF5

inputfile = "ModelParamIsochrone_roi.jl"

#CallAResponse.RunWmat(inputfile)

#CallAResponse.RunGfunc(inputfile)


include(inputfile)
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunM(inputfile,tabomega,VERBOSE=1)

#=

# compute the determinants with a gradient descent
bestomg = CallAResponse.FindZeroCrossing(inputfile,0.00,0.03,NITER=16,VERBOSE=1)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(inputfile,bestomg)
ModeR,ModeShape = CallAResponse.GetModeShape(inputfile,0.01,15.,100,EM)

=#
