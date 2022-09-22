
# need this to get the parameters...
inputfile = "ModelParamPlummer_ISO.jl"
include(inputfile)


import CallAResponse
using HDF5

# call the function to construct W matrices
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,
                      wmatdir,
                      FHT,
                      Kv,Kw,
                      basis,
                      lharmonic,
                      n1max,
                      Ω₀,
                      modelname,
                      rmin,rmax,
                      VERBOSE=2)



# call the function to compute G(u) functions
CallAResponse.RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,
                     ndFdJ,
                     wmatdir,gfuncdir,
                     FHT,
                     Kv,Kw,
                     basis,
                     lharmonic,
                     n1max,
                     Ω₀,
                     modelname,dfname,
                     rmin,rmax,
                     VERBOSE=1)


# construct a grid of frequencies to probe
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

# compute the matrix response at each location
tabdet = CallAResponse.RunM(tabomega,
                         dψ,d2ψ,
                         gfuncdir,modedir,
                         FHT,
                         Kv,Kw,
                         basis,
                         lharmonic,
                         n1max,
                         Ω₀,
                         modelname,dfname,
                         rmin,rmax,
                         VERBOSE=1)
