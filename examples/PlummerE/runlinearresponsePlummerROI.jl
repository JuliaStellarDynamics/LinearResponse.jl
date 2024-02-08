"""
for the radially-biased Plummer model: compute linear response theory
"""

# Bring in all parameters from the input file
inputfile = "ModelParamPlummerROIFiducial.jl"
inputfile = "ModelParamPlummerROIExample.jl"
include(inputfile)

import LinearResponse
using HDF5
using Plots

# package the Linear Response steps to compute M:
# 1. Call the function to construct W matrices
# 2. Run the G function calculation
# 3. Compute the matrix coefficients
@time LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
tabω = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
@time tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabω,FHT,Parameters)
@time tabdet = LinearResponse.RunDeterminant(tabω,FHT,Parameters)



tabOmega = collect(range(Omegamin,Omegamax,length=nOmega))
tabEta = collect(range(Etamin,Etamax,length=nEta))
    
epsilon = abs.(reshape(tabdet,nEta,nOmega))
println(minimum(epsilon))

# Plot
contour(tabOmega,tabEta,log10.(epsilon), levels=10, color=:black, #levels=[-2.0, -1.5, -1.0, -0.5, -0.25, 0.0], 
        xlabel="Re[ω]", ylabel="Im[ω]", xlims=(Omegamin,Omegamax), ylims=(Etamin,Etamax),
        clims=(-2, 0), aspect_ratio=:equal, legend=false)
savefig("ROIdeterminant.png")



# find a pole by using gradient descent
startingomg = 0.0 + 0.05im
@time bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
