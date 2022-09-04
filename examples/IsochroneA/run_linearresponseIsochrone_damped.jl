
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
                               Ω0,
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
                                Ω0,
                                modelname,dfname,
                                rb,
                                bc=bc,G=G,M=M,
                                VERBOSE=1)


#tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabdet = CallAResponse.RunMIsochrone(inputfile,tabomega,VERBOSE=0)

# for the minimum, go back and compute the mode shape
#bestomg = 0.0143 - 0.00141im # for n1max = 10
#bestomg = 0.01 - 0.000733im # for n1max = 12
bestomg = 0.00838706563046674 - 0.0005277615331419046im
EV,EF,EM = CallAResponse.ComputeModeTablesIsochrone(inputfile,bestomg)
ModeR,ModeShape = CallAResponse.GetModeShape(inputfile,0.01,25.,100,EM)
