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
G  = 1.
rb = 8.0
lmax,nmax = 2,16 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

#####
# Model Potential
#####


modelname = "IsochroneA"

bc, M, G = 1.,1.,1.
potential(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)



function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    return ISOndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG)
end



function ISOndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # Function that returns the value of n.dF/dJ where the DF follows the normalisation convention int dxdv F = Mtot

    dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG) # Current value of dF/E. ATTENTION, the DF is assumed to be isotropic
    res = ndotOmega*dFdE # Current value of n.dF/dJ. ATTENTION, the DF is assumed to be isotropic
    #####
    return res
end


#####
# Parameters
#####
K_u        = 150    # number of Legendre integration sample points
K_v        = 100    # number of allocations is directly proportional to this
NstepsWMat = 100    # number of allocations is insensitive to this (also time, largely?

lharmonic = 1
n1max = 10  # maximum number of radial resonances to consider

# Mode of response matrix computation
LINEAR = "damped"


#####
# Outputs directories
#####
wmatdir="wmat/"
gfuncdir="gfunc/"
modedir = "xifunc/"

# WARNING : / at the end to check !
