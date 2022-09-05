
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochrone_roi.jl"
include(inputfile)

import CallAResponse
using HDF5

# compute the Fourier-transformed basis elements
CallAResponse.RunWmatIsochrone(wmatdir,
                               K_u,K_v,K_w,
                               basis,
                               lharmonic,
                               n1max,
                               nradial,
                               Ω0,
                               modelname,
                               rb,
                               VERBOSE=2)

# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(ndFdJ,
                                wmatdir,gfuncdir,
                                K_u,K_v,K_w,
                                basis,
                                lharmonic,
                                n1max,
                                nradial,
                                Ω0,
                                modelname,dfname,
                                rb,
                                VERBOSE=1)

# make a grid of omegas to test
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabdet = CallAResponse.RunMIsochrone(tabomega,
                                     gfuncdir,modedir,
                                     K_u,K_v,K_w,
                                     basis,
                                     lharmonic,
                                     n1max,
                                     nradial,
                                     Ω0,
                                     modelname,dfname,
                                     rb,
                                     VERBOSE=1)

# compute the determinants with a gradient descent to find the exact frequency for ROI
bestomg = CallAResponse.FindZeroCrossingIsochrone(0.00,0.03,ψ,dψ,d2ψ,
                                         gfuncdir,modedir,
                                         K_u,K_v,K_w,
                                         basis,
                                         lharmonic,
                                         n1max,
                                         nradial,
                                         Ω0,
                                         modelname,dfname,
                                         rb,
                                         rmin,rmax,NITER=16,VERBOSE=1)

bestomg = 0.0 + 0.022im

# for the minimum, go back and compute the mode shape:
# need an isochrone specific version
EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,ψ,dψ,d2ψ,
                                          gfuncdir,modedir,
                                          K_u,K_v,K_w,
                                          basis,
                                          lharmonic,
                                          n1max,
                                          nradial,
                                          Ω0,
                                          modelname,dfname,
                                          rb,
                                          VERBOSE=1)

ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,lharmonic,
                                              0.01,15.,100,EM,VERBOSE=1)
