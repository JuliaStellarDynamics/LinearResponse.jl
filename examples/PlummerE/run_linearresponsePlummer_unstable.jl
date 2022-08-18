"""
for the Plummer model: compute linear response theory
"""

import CallAResponse
using HDF5

inputfile = "ModelParamPlummer_ROI.jl"

# compute the Fourier-transformed basis functions
#CallAResponse.RunWmat(inputfile)

# compute the G functions
CallAResponse.RunGfunc(inputfile)


include(inputfile)
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunM(inputfile,tabomega,VERBOSE=1)

#=
# compute the determinants with a gradient descent
bestomg = CallAResponse.FindZeroCrossing(inputfile,0.00,0.03,NITER=16,VERBOSE=1)

# for n1max=2
#bestomg =

println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(inputfile,bestomg)
ModeR,ModeShape = CallAResponse.GetModeShape(inputfile,0.01,15.,100,EM)

=#
