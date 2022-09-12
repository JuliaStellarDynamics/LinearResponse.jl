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
-modelname

"""


import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


##############################
# Basis
##############################
G  = 1.

# Clutton-Brock (1972) basis
basisname = "CluttonBrock"
rb = 14.
lmax,nmax = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB72Basis_create(lmax=lmax,nmax=nmax,G=G,rb=rb) 

# # Kalnajs (1976) basis
# basisname = "Kalnajs"
# rb, kKA = 5., 7
# lmax,nmax = 2,7
# basis = AstroBasis.K76Basis_create(lmax=lmax,nmax=nmax,G=G,rb=rb,kKA=kKA)

##############################
# Model Potential
##############################
const modelname = "Mestel"

const R0, V0 = 20., 1.
const ψ(r::Float64)::Float64   = OrbitalElements.ψMestel(r,R0,V0)
const dψ(r::Float64)::Float64  = OrbitalElements.dψMestel(r,R0,V0)
const d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψMestel(r,R0,V0)
const d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψMestel(r,R0,V0)
const d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψMestel(r,R0,V0)
const Ω₀ = OrbitalElements.Ω₀Mestel(R0,V0)
println("Ω₀ = ",Ω₀)

##############################
# Outputs directories
##############################
const wmatdir="wmat/"*basisname*"/"
const gfuncdir="gfunc/"*basisname*"/"
const modedir = "xifunc/"*basisname*"/"

##############################
# Model DF
##############################
const dfname = "Sellwood"
const q0 = 11.44
const σ0 = OrbitalElements.sigmar_Mestel_DF(R0,V0,q0)
const C0 = OrbitalElements.normC_Mestel_DF(G,R0,V0,q0)

const Rin, Rout, Rmax = 1., 11.5, 20. # Tapering radii
const xi=0.5                          # Self-gravity fraction
const nu, mu = 4, 5                   # Tapering exponants

const ndFdJ(n1::Int64,n2::Int64,
        E::Float64,L::Float64,
        ndotOmega::Float64)::Float64   = OrbitalElements.mestel_Zang_ndDFdJ(n1,n2,E,L,ndotOmega;
                                                                            R0=R0,Rin=Rin,Rmax=Rmax,
                                                                            V0=V0,
                                                                            xi=xi,C=C0,
                                                                            q=q0,sigma=σ0,
                                                                            nu=nu,mu=mu)


#####
# Parameters
#####
# Radii for frequency truncations
const rmin = 0.15
const rmax = 13.5

const Ku = 100           # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

const Kv = 100    # number of allocations is directly proportional to this
const Kw = 200    # number of allocations is insensitive to this (also time, largely?

const lharmonic = 2
const n1max = 4  # maximum number of radial resonances to consider