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
import PerturbPlasma
using HDF5


#####
# Basis
#####
G         = 1.     # the gravitational constant
rb        = 20.0   # the scale for the basis elements
lmax,nmax = 2,100  # usually lmax corresponds to the considered harmonics lharmonic
basis     = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim      = basis.dimension
nradial   = basis.nmax

println("Using dimension $ndim.")

#####
# Model Potential
#####
modelname = "IsochroneA"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Omega0 = OrbitalElements.Omega0Isochrone(bc,M,G)


#####
# Distribution Function
#####
dfname = "isotropic"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
               """caution, this is ISOTROPIC specific."""

   # Current value of dF/E
   dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG)

   # Current value of n.dF/dJ
   res = ndotOmega*dFdE

   return res
end


#####
# Outputs directories
#####
wmatdir    = "wmat/"
gfuncdir   = "gfunc/"
modedir    = "xifunc/"


#####
# Parameters
#####
K_u = 200      # number of Legendre integration sample points
K_v = 200      # number of allocations is directly proportional to this
K_w = 200      # number of allocations is insensitive to this (also time, largely?)

lharmonic  = 1        # which harmonic we are considering
n1max      = 10       # maximum number of radial resonances to consider


# Frequencies to probe
LINEAR     = "damped"
nOmega = 51
Omegamin = 0.0
Omegamax = 0.05
nEta = 50
Etamin = -0.005
Etamax = 0.0

#= # to see the computation in the upper half plane...
# Frequencies to probe
LINEAR     = "unstable"
nOmega = 61
Omegamin = -0.05
Omegamax = 0.05
nEta = 30
Etamin = 0.0001
Etamax = 0.005
=#
