
# need this to get the parameters...
inputfile = "ModelParamPlummer_ISO.jl"
include(inputfile)

import CallAResponse
using HDF5

# call the function to construct W matrices
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,
                      wmatdir,
                      K_u,K_v,K_w,
                      basis,
                      lharmonic,
                      n1max,
                      nradial,
                      Ω0,
                      modelname,
                      rb,
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
                       Ω0,
                       modelname,dfname,
                       rb,
                       rmin,rmax,
                       VERBOSE=1)


# construct a grid of frequencies to probe
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

# compute the matrix response at each location
tabdet = CallAResponse.RunM(tabomega,
                            dψ,d2ψ,
                            gfuncdir,modedir,
                            K_u,K_v,K_w,
                            basis,
                            lharmonic,
                            n1max,
                            nradial,
                            Ω0,
                            modelname,dfname,
                            rb,
                            rmin,rmax,
                            VERBOSE=1)


# compute the determinants with a gradient descent
bestomg = CallAResponse.FindZeroCrossing(0.00,0.03,dψ,d2ψ,
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

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.ComputeModeTables(bestomg,dψ,d2ψ,
                                         gfuncdir,modedir,
                                         K_u,K_v,K_w,
                                         basis,
                                         lharmonic,
                                         n1max,
                                         nradial,
                                         Ω0,
                                         modelname,dfname,
                                         rb,
                                         rmin,rmax,VERBOSE=1)

ModeR,ModeShape = CallAResponse.GetModeShape(basis,lharmonic,
                                             0.01,15.,100,EM,VERBOSE=1)
