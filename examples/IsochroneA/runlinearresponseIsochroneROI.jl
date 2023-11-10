
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneROI.jl"
include(inputfile)

import LinearResponse
using HDF5

println(Parameters.nbResVec)
println(Parameters.tabResVec)

# hack to only compute with one resonance vector
Parameters.nbResVec = 1
Parameters.tabResVec = zeros(Int64,2,Parameters.nbResVec)
Parameters.tabResVec[1,1], Parameters.tabResVec[2,1] = -1, 2

println(Parameters.nbResVec)
println(Parameters.tabResVec)

Parameters.Ku = 250
Parameters.Kv = 350
Parameters.Kw = 200

FHT = FiniteHilbertTransform.LegendreFHT(Parameters.Ku)
#@time LinearResponse.RunWmatIsochrone(FHT,bc,M,G,basis,Parameters)
@time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

"""
# compute the Fourier-transformed basis elements
#LinearResponse.RunWmatIsochrone(FHT,bc,M,G,basis,Parameters)

# compute the G(u) functions
LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# compute the matrix response at each location in tabomega
tabdet = LinearResponse.RunAXi(FHT,Parameters)

# construct a grid of frequencies to probe
#tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.02im
#startingomg = 0.0 - 0.0000011773114542602149im
bestomg = LinearResponse.FindPole(startingomg,FHT,Parameters,1.e-12)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 200
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
"""
