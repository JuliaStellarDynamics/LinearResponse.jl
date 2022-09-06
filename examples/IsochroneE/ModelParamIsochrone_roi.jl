"""
an example input file for running all steps in estimating the Linear Response for a given model

Must include:
-ψ
-dψ
-d2ψ
-basis
-ndim
-nradial
-ndFdJ


"""


import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
using HDF5


#####
# Basis
#####
G  = 1.
#rb = 20.0
rb = 15.0
lmax,nmax = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

# Model Potential
modelname = "IsochroneE"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)

rmin = 0.0
rmax = 10000.0


dfname = "roi1.0"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)

    Q = OrbitalElements.isochrone_Q_ROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochrone_Saha_dDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochrone_dQdE_ROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochrone_dQdL_ROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end



# integration parameters
Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

lharmonic = 2
n1max = 2  # maximum number of radial resonances to consider

# Mode of response matrix computation
# Frequencies to probe
nOmega = 51
Omegamin = -0.02
Omegamax = 0.02
nEta = 50
Etamin = 0.001
Etamax = 0.04



# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"



# WARNING : / at the end to check !
