"""
for the radially-biased Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerROIhires.jl"
inputfile = "ModelParamPlummerROIFiducial.jl"

include(inputfile)

import LinearResponse
using HDF5

"""
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
@time LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)
@time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

"""

# call the function to construct W matrices
LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)

# run the G function calculation
LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# compute the matrix response at each location in tabomega
tabdet = LinearResponse.RunAXi(FHT,Parameters)

# combine all three steps into one wrapper:
#LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
#tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.1 + 0.01im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)#,1.e-12)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
