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
