"""
for the Isochrone radially-unstable model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamIsochroneROIhires.jl"
include(inputfile)

import LinearResponse
using HDF5

println(Parameters.nbResVec)
println(Parameters.tabResVec)

"""
# hack to only compute with one resonance vector
Parameters.nbResVec = 1
Parameters.tabResVec = zeros(Int64,2,Parameters.nbResVec)
Parameters.tabResVec[1,1], Parameters.tabResVec[2,1] = -1, 2

println(Parameters.nbResVec)
println(Parameters.tabResVec)

Kutest = [25,50,100,200,300,400,500,600,700,800,1600,3200]
Kvtest = [100,200,300,400,500,600,700,800]
#Kutest = [101,201,301,401,501,601,701,801] # test with VMAPN = 1
"""
#Parameters.Ku = 200
#Parameters.Kv = 200
#Parameters.Kw = 200

#FHT = FiniteHilbertTransform.LegendreFHT(Parameters.Ku)
#@time LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)
#@time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

"""
for Kval in Kutest
    #Parameters.Ku = Kval
    #Parameters.Ku = 400 # for Kv tests
    #Parameters.Kv = Kval
    FHT = FiniteHilbertTransform.LegendreFHT(Parameters.Ku)
    @time LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)
    @time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)
end
"""



# combine all steps into one wrapper:
LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.01im
bestomg = LinearResponse.FindPole(startingomg,FHT,Parameters)#,1.e-12)

bestomg = 0.0 + 0.023132620586365966im#0.02271406012170436im
println("The zero-crossing frequency is #bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
