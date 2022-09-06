
# this input file is set up to mimic the Fouvry & Prunet (2022) damped mode result exactly
inputfile = "ModelParamIsochrone_damped.jl"

# need this to get the parameters...
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
                               Ω₀,
                               modelname,
                               rb,
                               bc=bc,G=G,M=M
                               VERBOSE=2)

# compute the G(u) functions
CallAResponse.RunGfuncIsochrone(ndFdJ,
                                wmatdir,gfuncdir,
                                K_u,K_v,K_w,
                                basis,
                                lharmonic,
                                n1max,
                                nradial,
                                Ω₀,
                                modelname,dfname,
                                rb,
                                bc=bc,G=G,M=M,
                                VERBOSE=1)


tabdet = CallAResponse.RunMIsochrone(tabomega,
                                     gfuncdir,modedir,
                                     K_u,K_v,K_w,
                                     basis,
                                     FHT,
                                     lharmonic,
                                     n1max,
                                     nradial,
                                     Ω₀,
                                     modelname,dfname,
                                     rb,
                                     rmin,rmax,
                                     VERBOSE=1)

 # for the minimum, go back and compute the mode shape
 #bestomg = 0.0143 - 0.00141im # for n1max = 10
 #bestomg = 0.01 - 0.000733im # for n1max = 12
 bestomg = 0.00838706563046674 - 0.0005277615331419046im


 # for the minimum, go back and compute the mode shape
 EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,ψ,dψ,d2ψ,
                                            gfuncdir,modedir,
                                            K_u,K_v,K_w,
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
