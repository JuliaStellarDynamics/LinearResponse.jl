"""
for the Plummer model: compute linear response theory
"""


# bring in all parameters from the input file
inputfile = "ModelParamPlummer_ROI.jl"
include(inputfile)

import CallAResponse
using HDF5

#


#Kv = 1600
n1max = 1
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



for aval=1:7

   raval = (aval-1)*0.05 + 0.75
   dfnamein = "roi"*string(raval)

   #println("running on $dfnamein")

   ndFdJin(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64) = ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc=bc,M=M,astronomicalG=G,Ra=raval)

   # call the function to compute G(u) functions
   CallAResponse.RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,
                          ndFdJin,
                          wmatdir,gfuncdir,
                          FHT,
                          Kv,Kw,
                          basis,
                          lharmonic,
                          n1max,
                          Ω₀,
                          modelname,dfnamein,
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
                                modelname,dfnamein,
                                rmin,rmax,
                                VERBOSE=-1)

    startingomg = 0.0 + 0.2im;
    bestomg = CallAResponse.FindPole(startingomg,FHT,lharmonic,n1max,nradial,ndim,dψ,d2ψ,Ω₀,rmin,rmax,modedir,modelname,dfnamein,rb,Kv,VERBOSE=4,KuTruncation=10)


end
