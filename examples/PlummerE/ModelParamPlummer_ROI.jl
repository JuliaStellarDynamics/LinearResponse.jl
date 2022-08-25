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
import PerturbPlasma
using HDF5


#####
# Basis
#####
G  = 1.
rb = 20.0
lmax,nmax = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

#####
# Model Potential
#####


const modelname = "PlummerE"

const bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.plummer_dddpsi_dddr(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddddpsi_ddddr(r,bc,M,G)
Omega0 = OrbitalElements.plummer_Omega0(bc,M,G)


#dfname = "roi1.0"
#dfname = "roi0.75"
dfname = "roi0.90"
#dfname = "roi1.1"


function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=0.9)

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# integration parameters
K_u = 200    # number of Legendre integration sample points
K_v = 200    # number of allocations is directly proportional to this
K_w = 200    # number of allocations is insensitive to this (also time, largely?

lharmonic = 2
n1max     = 2 # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
xifuncdir= "xifunc/"
modedir  = "xifunc/"

#=
# Frequencies to probe
LINEAR   = "unstable"
nOmega   = 51
Omegamin = -0.02
Omegamax = 0.02
nEta     = 50
Etamin   = 0.001
Etamax   = 0.05
=#

# Frequencies to probe
LINEAR   = "damped"
nOmega   = 51
Omegamin = -0.02
Omegamax = 0.02
nEta     = 50
Etamin   = -0.005
Etamax   = -0.0001




# WARNING : / at the end to check !
