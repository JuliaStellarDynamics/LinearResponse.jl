"""
an input file which reproduces the Fouvry & Prunet (2022) unstable l=2 mode calculation
driven by run_linearresponseIsochrone_unstable.jl
"""

import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


# choose a basis for computation of the Fourier-transformed basis elements
G  = 1.
rb = 20.0
lmax,nmax = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
basis     = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim      = basis.dimension
nradial   = basis.nmax

# choose a model potential
modelname = "IsochroneA"

bc, M, G = 1.,1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Omega0 = OrbitalElements.Omega0Isochrone(bc,M,G)


# choose a distribution function for G(u) calculation
# if you change this, you will have to re-run the G calculation!
dfname = "roi1.0"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)

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


# integration arameters
K_u = 200    # number of Legendre integration sample points
K_v = 200    # number of allocations is directly proportional to this
K_w = 200    # number of allocations is insensitive to this (also time, largely?

lharmonic = 2
n1max     = 2  # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

# Frequencies to probe
LINEAR   = "unstable"
nOmega   = 61
Omegamin = -0.05
Omegamax = 0.05
nEta     = 30
Etamin   = 0.0001
Etamax   = 0.005


# WARNING : / at the end to check !
