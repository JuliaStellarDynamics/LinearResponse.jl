"""
for the isotropic Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerISO.jl"
inputfile = "ModelParamPlummerISOFiducial.jl"

include(inputfile)

import LinearResponse
using HDF5

# combine all three steps into one wrapper:
LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# call the function to construct W matrices
#LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)

# run the G function calculation
#LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# compute the matrix response at each location in tabomega
#tabdet = LinearResponse.RunAXi(FHT,Parameters)


# construct a grid of frequencies to probe
#tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
#tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)
#tabdetXi = LinearResponse.RunDeterminant(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.45 + 0.002im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)

#bestomg = 0.0 + 0.02271406012170436im
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)
println(EM) # print the eigenvalues just to check

modeRmin = 0.01
modeRmax = 10.0
nmode = 400
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
