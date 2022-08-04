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
lmax,nmax = 1,100  # usually lmax corresponds to the considered harmonics lharmonic
basis     = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim      = basis.dimension
nradial   = basis.nmax

#####
# Model Potential
#####
modelname = "IsochroneA"

bc, M = 1.,1.
potential(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)


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
# Parameters
#####
K_u        = 200      # number of Legendre integration sample points
K_v        = 200      # number of allocations is directly proportional to this
NstepsWMat = 200      # number of allocations is insensitive to this (also time, largely?)

lharmonic  = 1        # which harmonic we are considering
n1max      = 10       # maximum number of radial resonances to consider

LINEAR     = "damped" # Mode of response matrix computation

#####
# Outputs directories
#####
wmatdir    = "wmat/"
gfuncdir   = "gfunc/"
modedir    = "xifunc/"
