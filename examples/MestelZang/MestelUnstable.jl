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
import FiniteHilbertTransform
using HDF5


##############################
# Basis
##############################
const G  = 1.

# Clutton-Brock (1972) basis
const basisname = "CluttonBrock"
const rb = 20.0
const lmax,nmax = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
const basis = AstroBasis.CB72Basis_create(lmax=lmax,nmax=nmax,G=G,rb=rb) 

# # Kalnajs (1976) basis
# basisname = "Kalnajs"
# rb, kKA = 5., 7
# lmax,nmax = 2,2
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
const q0 = 6
const σ0 = OrbitalElements.sigmar_Mestel_DF(R0,V0,q0)
const C0 = OrbitalElements.normC_Mestel_DF(G,R0,V0,q0)

const Rin, Rout, Rmax = 1., 11.5, 20. # Tapering radii
const xi=1.0                         # Self-gravity fraction
const nu, mu = 6, 5                   # Tapering exponants

const dfname = "Zang_q_"*string(q0)*"_xi_"*string(xi)*"_mu_"*string(mu)*"_nu_"*string(nu)

const ndFdJ(n1::Int64,n2::Int64,
        E::Float64,L::Float64,
        ndotOmega::Float64)::Float64   = OrbitalElements.mestel_Zang_ndDFdJ(n1,n2,E,L,ndotOmega;
                                                                            R0=R0,Rin=Rin,Rmax=Rmax,
                                                                            V0=V0,
                                                                            xi=xi,C=C0,
                                                                            q=q0,sigma=σ0,
                                                                            nu=nu,mu=mu)

##############################
# Legendre / Integration parameters
##############################
# Radii for frequency truncations
const rmin = 0.1
const rmax = 20.0

const Ku = 100           # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

const Kv = 100    # number of allocations is directly proportional to this
const Kw = 200    # number of allocations is insensitive to this (also time, largely?

const lharmonic = 2
const n1max = 1  # maximum number of radial resonances to consider