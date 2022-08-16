


import CallAResponse
using HDF5

# to mimic Fouvry & Prunet exactly
#inputfile = "ModelParamIsochrone_roi.jl"

# create a reduced version for checking
inputfile = "ModelParamIsochrone_roi_reduced.jl"

# compute the Fourier-transformed basis elements
CallAResponse.RunWmatIsochrone(inputfile)

# compute the G(u) functions
#CallAResponse.RunGfuncIsochrone(inputfile)

# compute the determinants with a gradient descent
#bestomg = CallAResponse.FindZeroCrossing(inputfile,0.00,0.03,NITER=16)
#println(bestomg)

# for the minimum, go back and compute the mode shape
#EV,EF,EM = CallAResponse.RunShape(inputfile,bestomg)
#ModeR,ModeShape = CallAResponse.tabShapeMode(inputfile,0.01,15.,100,EM)
