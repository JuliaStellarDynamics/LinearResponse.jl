"""

all steps combined into one: could pause and restart if this was too much.


"""

import LinearResponse

inputfile = "SellwoodStable.jl"
include(inputfile)

# call the function to construct W matrices
LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,params)

# call the function to compute G(u) functions
LinearResponse.RunGfunc(ndFdJ,FHT,params)

# call the function to compute Xi decomposition coefficients
LinearResponse.RunAXi(FHT,params)

# construct a grid of frequencies to probe
nbω0 = 100                 # Number of ω0 for which the matrix is computed
ω0min, ω0max = 0., 1.5     # Minimum and maximum ω0
nbη = 50                  # Number of η for which the matrix is computed
ηmin, ηmax = -0.1, 0.3      # Minimum and maximum η
tabω = LinearResponse.gridomega(ω0min,ω0max,nbω0,ηmin,ηmax,nbη)
# compute the matrix response at each location
tabdet = LinearResponse.RunDeterminant(tabω,FHT,params)

# Mode Finding
Ωguess = 1.0
ηguess = 0.0
ωguess = Ωguess + im*ηguess
ωMode = LinearResponse.FindPole(ωguess,FHT,params)
println("ωMode = ",ωMode)

# # Mode Shape
# ωMode = 0.28 - 0.017*im
# EV,EF,EM = CallAResponse.ComputeModeTables(ωMode,FHT,basis,params)
# ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShape(basis,0.1,10.,1000,EF,params)