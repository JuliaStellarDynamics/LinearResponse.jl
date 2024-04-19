"""

all steps combined into one: could pause and restart if this was too much.


"""

using LinearResponse
using Plots

inputfile = "SellwoodStable.jl"
inputfile = "MestelUnstable.jl"

include(inputfile)

# call the function to construct W matrices
LinearResponse.RunWmat(model,FHT,basis,params)

# call the function to compute G(u) functions
LinearResponse.RunGfunc(distributionfunction,FHT,params)

# call the function to compute Xi decomposition coefficients
LinearResponse.RunAXi(FHT,params)

# # construct a grid of frequencies to probe
nbω0 = 100                 # Number of ω0 for which the matrix is computed
ω0min, ω0max = 0., 2.5     # Minimum and maximum ω0
nbη = 50                  # Number of η for which the matrix is computed
ηmin, ηmax = -0.1, 0.6      # Minimum and maximum η
# tabω = LinearResponse.gridomega(ω0min,ω0max,nbω0,ηmin,ηmax,nbη)
# # compute the matrix response at each location
# tabdet = LinearResponse.RunDeterminant(tabω,FHT,params)

# construct a grid of frequencies to probe
tabω = LinearResponse.gridomega(ω0min,ω0max,nbω0,ηmin,ηmax,nbη)
@time tabdet = LinearResponse.RunDeterminant(tabω,FHT,params)


tabOmega = collect(range(ω0min,ω0max,length=nbω0))
tabEta = collect(range(ηmin,ηmax,length=nbη))
    
epsilon = abs.(reshape(tabdet,nbη,nbω0))

# Plot
contour(tabOmega,tabEta,log10.(epsilon), levels=10, color=:black, #levels=[-2.0, -1.5, -1.0, -0.5, -0.25, 0.0], 
        xlabel="Re[ω]", ylabel="Im[ω]", xlims=(ω0min,ω0max), ylims=(ηmin,ηmax),
        clims=(-2, 0), aspect_ratio=:equal, legend=false)
savefig("ROIdeterminant.png")



# Mode Finding
Ωguess = 1.0
ηguess = 0.5
ωguess = Ωguess + im*ηguess
ωMode = LinearResponse.FindPole(ωguess,FHT,params)
println("ωMode = ",ωMode)

# # Mode Shape
# ωMode = 0.28 - 0.017*im
# EV,EF,EM = CallAResponse.ComputeModeTables(ωMode,FHT,basis,params)
# ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,0.1,10.,1000,EF,params)