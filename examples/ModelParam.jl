import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


#####
# Basis
#####
G  = 1.
rb = 1.
lmax,nmax = 2,10 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB72Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

#####
# Model Potential
#####

#=
modelname = "isochroneE"

bc, M, G = 1.,1.,1.
potential(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
=#
#=
modelname = "PlummerE"

bc, M, G = 1.,1.,1.
potential(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.plummer_Omega0(bc,M,G)
=#

modelname = "MestelZ"

R0, V0 = 20., 1.
eps0 = 0.001
potential(r::Float64)::Float64   = OrbitalElements.mestel_psi(r,R0,V0,eps0)
dpotential(r::Float64)::Float64  = OrbitalElements.mestel_dpsi_dr(r,R0,V0,eps0)
ddpotential(r::Float64)::Float64 = OrbitalElements.mestel_ddpsi_ddr(r,R0,V0,eps0)
Omega0 = OrbitalElements.mestel_Omega0(R0,V0,eps0)

#####
# Model DF
#####
q0 = 11.44
sigma0 = OrbitalElements.sigmar_Mestel_DF(R0,V0,q0)
C0 = OrbitalElements.normC_Mestel_DF(R0,V0,q0)

Rin, Rout, Rmax = 1., 11.5, 20. # Tapering radii
xi=0.5                          # Self-gravity fraction
nu, mu = 4, 5                   # Tapering exponants

ndFdJ(n1::Int64,n2::Int64,
        E::Float64,L::Float64,
        ndotOmega::Float64)::Float64   = OrbitalElements.mestel_Zang_ndDFdJ(n1,n2,E,L,ndotOmega;
                                                                            R0=R0,Rin=Rin,Rmax=Rmax,
                                                                            V0=V0,
                                                                            xi=xi,C=C0,
                                                                            q=q0,sigma=sigma0,
                                                                            nu=nu,mu=mu)

#####
# Parameters
#####
K_u = 150           # number of Legendre integration sample points
K_v        = 100    # number of allocations is directly proportional to this
NstepsWMat = 100    # number of allocations is insensitive to this (also time, largely?

lharmonic = 2 
n1max = 4  # maximum number of radial resonances to consider

# Mode of response matrix computation
LINEAR = "damped"


#####
# Outputs directories
#####
wmatdir="wmat/"
gfuncdir="gfunc/"

# WARNING : / at the end to check !