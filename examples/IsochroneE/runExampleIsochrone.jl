"""
for the radially-biased Plummer model: compute linear response theory
"""

using OrbitalElements
using AstroBasis
using FiniteHilbertTransform
using LinearResponse
using HDF5

using Plots


# Basis
G  = 1.
rb = 2.0
lmax,nradial = 2,5 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)



# Model Potential
const modelname = "IsochroneE2"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.NumericalIsochrone()

rmin = 0.0
rmax = Inf


dfname = "roi1.0"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.0)

    Q = OrbitalElements.isochroneQROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochroneSahadDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochronedQdEROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochronedQdLROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end


# Linear Response integration parameters
Ku = 50    # number of Legendre integration sample points
Kv = 30    # number of allocations is directly proportional to this
Kw = 20    # number of allocations is insensitive to this (also time, largely)?


# Define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max = 1 

# output directories
wmatdir  = "./"
gfuncdir = "./"
modedir  = "./"

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 280
Omegamin = -0.2
Omegamax = 0.2
nEta     = 280
Etamin   = -0.1
Etamax   = 0.1


VERBOSE   = 2
OVERWRITE = true
VMAPN     = 1
ADAPTIVEKW= false

OEparams = OrbitalElements.OrbitalParameters(EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ω₀=OrbitalElements.frequency_scale(model),Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)



# package the Linear Response steps to compute M:
# 1. Call the function to construct W matrices
# 2. Run the G function calculation
# 3. Compute the matrix coefficients
#@time LinearResponse.RunLinearResponse(model,ndFdJ,FHT,basis,Parameters)

@time LinearResponse.RunWmat(model,FHT,basis,Parameters)

# call the function to compute G(u) functions
@time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# call the function to compute decomposition coefficients
@time LinearResponse.RunAXi(FHT,Parameters)

#MMat, tabaMcoef, tabωminωmax = LinearResponse.PrepareM(Parameters)

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
        clims=(-2, 0), aspect_ratio=:equal, legend=false)
savefig("ROIdeterminant2.png")


# find a pole by using gradient descent
startingomg = 0.0 + 0.01im
@time bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

if isfinite(bestomg)
        # for the minimum, go back and compute the mode shape
        EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

        modeRmin = 0.01
        modeRmax = 15.0
        nmode = 100
        ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
else
        println("runExampleIsochrone.jl: no mode found.")
end