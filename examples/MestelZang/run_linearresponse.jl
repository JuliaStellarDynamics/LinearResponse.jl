"""

all steps combined into one: could pause and restart if this was too much.


"""

import CallAResponse

inputfile = "MestelUnstable.jl"
include(inputfile)

const verbose = 1
const overwrite = false
# call the function to construct W matrices
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,wmatdir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,rmin,rmax,VERBOSE=verbose)

# # call the function to compute G(u) functions
CallAResponse.RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,wmatdir,gfuncdir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,dfname,rmin,rmax,VERBOSE=verbose,OVERWRITE=overwrite)

# construct a grid of frequencies to probe
nbω0 = 50                 # Number of ω0 for which the matrix is computed
ω0min, ω0max = 0.5, 1.5     # Minimum and maximum ω0
nbη = 50                   # Number of ω0 for which the matrix is computed
ηmin, ηmax = -0.1, 0.5      # Minimum and maximum ω0
tabomega = CallAResponse.gridomega(ω0min,ω0max,nbω0,ηmin,ηmax,nbη)
# compute the matrix response at each location
tabdet = CallAResponse.RunM(tabomega,dψ,d2ψ,gfuncdir,modedir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,dfname,rmin,rmax,VERBOSE=1,OVERWRITE=overwrite)

# Mode Finding
Ωguess = 0.9
ηguess = 0.2
ωguess = Ωguess + im*ηguess
nradial, ndim, rb = basis.nmax, basis.dimension, basis.rb
ωMode = CallAResponse.FindPole(ωguess,FHT,lharmonic,n1max,nradial,ndim,dψ,d2ψ,Ω₀,rmin,rmax,modedir,modelname,dfname,rb,VERBOSE=verbose)
println("ωMode = ",ωMode)

# Mode Shape
EV,EF,EM = CallAResponse.ComputeModeTables(ωMode,dψ,d2ψ,FHT,gfuncdir,modedir,Kv,basis,lharmonic,n1max,Ω₀,rmin,rmax,modelname,dfname,VERBOSE=verbose)
ModeRadius,ModePotentialShape,ModeDensityShape = CallAResponse.GetModeShapeComplex(basis,lharmonic,n1max,0.1,20.,1000,EF,modedir,modelname,dfname,Ku,VERBOSE=1)