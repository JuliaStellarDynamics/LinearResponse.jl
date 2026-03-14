"""
for the isotropic Plummer model: compute linear response theory
"""

# bring in all parameters from the input file
inputfile = "ModelParamPlummerISOFiducial.jl"
include(inputfile)

using LinearResponse
using HDF5

LinearResponse.RunLinearResponse(model,distributionfunction,FHT,basis,Parameters)

# construct a grid of frequencies to probe
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabomega,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.05 + 0.01im
bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

# for the minimum, go back and compute the mode shape
EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)


modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)



MMat, tabaMcoef, tabωminωmax = LinearResponse.PrepareM(Parameters)

nOmega   = 280
Omegamin = -0.2
Omegamax = 0.2
nEta     = 280
Etamin   = -0.1
Etamax   = 0.1

# construct a grid of frequencies to probe
tabω = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)
@time tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabω,FHT,Parameters)
@time tabdet = LinearResponse.RunDeterminant(tabω,FHT,Parameters)


tabOmega = collect(range(Omegamin,Omegamax,length=nOmega))
tabEta = collect(range(Etamin,Etamax,length=nEta))
    
epsilon = abs.(reshape(tabdet,nEta,nOmega))

# Plot
contour(tabOmega,tabEta,log10.(epsilon), levels=10, color=:black, #levels=[-2.0, -1.5, -1.0, -0.5, -0.25, 0.0], 
        xlabel="Re[ω]", ylabel="Im[ω]", xlims=(Omegamin,Omegamax), ylims=(Etamin,Etamax),
        clims=(-5, 0), aspect_ratio=:equal, legend=false)
savefig("ISOdeterminant.png")


println("The zero-crossing frequency is $bestomg.")