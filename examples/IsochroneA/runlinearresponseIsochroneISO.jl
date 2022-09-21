
# this input file is set up to mimic the Fouvry & Prunet (2022) damped mode result exactly
inputfile = "ModelParamIsochroneISO.jl"

# need this to get the parameters...
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
                               VERBOSE=-2)



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
                                VERBOSE=-1)



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

 # for the minimum, go back and compute the mode shape
 #bestomg = 0.0143 - 0.00141im # for n1max = 10
 #bestomg = 0.01 - 0.000733im # for n1max = 12
 bestomg = 0.00838706563046674 - 0.0005277615331419046im


 # for the minimum, go back and compute the mode shape
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
