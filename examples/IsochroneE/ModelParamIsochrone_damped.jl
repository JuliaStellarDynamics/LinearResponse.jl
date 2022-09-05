"""
an example input file for running all steps in estimating the Linear Response for a given model

Must include:
-potential
-dpotential
-ddpotential
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
rb = 20.0
lmax,nmax = 1,10 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

#####
# Model Potential
#####


modelname = "IsochroneE"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω0 = OrbitalElements.Ω₀Isochrone(bc,M,G)

rmin = 0.0
rmax = 10000.0



dfname = "isotropic"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG) # Current value of dF/E. ATTENTION, the DF is assumed to be isotropic
    res = ndotOmega*dFdE # Current value of n.dF/dJ. ATTENTION, the DF is assumed to be isotropic
    #####
    return res
end


# integration parameters
K_u        = 200    # number of Legendre integration sample points
K_v        = 200    # number of allocations is directly proportional to this
K_w        = 200    # number of allocations is insensitive to this (also time, largely?


lharmonic = 1
n1max     = 10  # maximum number of radial resonances to consider

# output directories
wmatdir   = "wmat/"
gfuncdir  = "gfunc/"
modedir   = "xifunc/"

# Frequencies to probe
nOmega   = 51
Omegamin = 0.0
Omegamax = 0.05
nEta     = 50
Etamin   = -0.005
Etamax   = 0.0


# WARNING : / at the end to check !
