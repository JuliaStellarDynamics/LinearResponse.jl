


import CallAResponse
using HDF5

inputfile = "ModelParamIsochrone_damped.jl"

# compute the Fourier-transformed basis elements
#CallAResponse.RunWmatIsochrone(inputfile)

# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(inputfile)

# need to organise the omegalist here

# need this to get the parameters... is this the best fix?
include(inputfile)


tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunM(inputfile,tabomega,VERBOSE=1)


#=


det_filename(modedir,modelname,dfname,lharmonic,n1max,K_u)
=#
