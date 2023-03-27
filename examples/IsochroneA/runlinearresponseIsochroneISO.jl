
# this input file is set up to mimic the Fouvry & Prunet (2022) damped mode result exactly
inputfile = "ModelParamIsochroneISO.jl"

# need this to get the parameters...
include(inputfile)


import LinearResponse
using HDF5


"""


# compute the Fourier-transformed basis elements
LinearResponse.RunWmatIsochrone(wmatdir,
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
LinearResponse.RunGfuncIsochrone(ndFdJ,
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
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

tabdet = LinearResponse.RunMIsochrone(tabomega,
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
                                     VERBOSE=1,KuTruncation=30)

 # for the minimum, go back and compute the mode shape
 #bestomg = 0.0143 - 0.00141im # for n1max = 10
 #bestomg = 0.01 - 0.000733im # for n1max = 12
 bestomg = 0.00838706563046674 - 0.0005277615331419046im


 # for the minimum, go back and compute the mode shape
 EV,EF,EM = LinearResponse.ComputeModeTables(bestomg,ψ,dψ,d2ψ,
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

 ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,lharmonic,
                                                                             0.01,15.,100,EM,VERBOSE=1)

"""
