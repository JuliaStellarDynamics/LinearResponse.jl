
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
                       K_u,K_v,K_w,
                       basis,
                       lharmonic,
                       n1max,
                       nradial,
                       Ω₀,
                       modelname,dfname,
                       rb,
                       rmin,rmax,
                       VERBOSE=1)
