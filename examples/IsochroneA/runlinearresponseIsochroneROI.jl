
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneroi.jl"
include(inputfile)

import CallAResponse
using HDF5

# compute the Fourier-transformed basis elements
CallAResponse.RunWmatIsochrone(wmatdir,
                               Ku,Kv,Kw,
                               basis,
                               lharmonic,
                               n1max,
                               nradial,
                               Ω₀,
                               modelname,
                               rb,
                               VERBOSE=2)



# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(ndFdJ,
                                wmatdir,gfuncdir,
                                Ku,Kv,Kw,
                                basis,
                                lharmonic,
                                n1max,
                                nradial,
                                Ω₀,
                                modelname,dfname,
                                rb,
                                VERBOSE=1)


# make a grid of omegas to test
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunMIsochrone(tabomega,
                                     gfuncdir,modedir,
                                     Ku,Kv,Kw,
                                     basis,
                                     FHT,
                                     lharmonic,
                                     n1max,
                                     nradial,
                                     Ω₀,
                                     modelname,dfname,
                                     rb,
                                     VERBOSE=1)


startingomg = 0.0 + 0.02im;
bestomg = CallAResponse.FindPole(startingomg,FHT,lharmonic,n1max,nradial,ndim,dψ,d2ψ,Ω₀,rmin,rmax,modedir,modelname,dfname,rb,VERBOSE=2)


# for the minimum, go back and compute the mode shape:
# need an isochrone specific version
EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,ψ,dψ,d2ψ,
                                          gfuncdir,modedir,
                                          Ku,Kv,Kw,
                                          basis,
                                          lharmonic,
                                          n1max,
                                          nradial,
                                          Ω₀,
                                          modelname,dfname,
                                          rb,
                                          VERBOSE=1)

ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,lharmonic,
                                              0.01,15.,100,EM,VERBOSE=1)
