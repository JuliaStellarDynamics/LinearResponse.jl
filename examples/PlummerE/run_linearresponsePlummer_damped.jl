

import CallAResponse
using HDF5

inputfile = "ModelParamPlummer_ISO.jl"

# compute the Fourier-transformed basis elements
#CallAResponse.RunWmat(inputfile)

# compute the G(u) functions
#CallAResponse.RunGfunc(inputfile)

# need this to get the parameters...
include(inputfile)

tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunM(inputfile,tabomega,VERBOSE=1)
