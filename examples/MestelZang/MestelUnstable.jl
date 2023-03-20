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
import CallAResponse
using HDF5


##############################
# Basis
##############################
const G  = 1.

# Clutton-Brock (1972) basis
const basisname = "CluttonBrock"
const rb = 10.0
const lmax,nmax = 2,20 # Usually lmax corresponds to the considered harmonics lharmonic
const basis = AstroBasis.CB72BasisCreate(lmax=lmax,nmax=nmax,G=G,rb=rb) 

##############################
# Model Potential
##############################
const modelname = "Mestel"

const R0, V0 = 20., 1.
const ψ(r::Float64)   = OrbitalElements.ψMestel(r,R0,V0)
const dψ(r::Float64)  = OrbitalElements.dψMestel(r,R0,V0)
const d2ψ(r::Float64) = OrbitalElements.d2ψMestel(r,R0,V0)
const d3ψ(r::Float64) = OrbitalElements.d3ψMestel(r,R0,V0)
const d4ψ(r::Float64) = OrbitalElements.d4ψMestel(r,R0,V0)
const Ω₀ = OrbitalElements.Ω₀Mestel(R0,V0)

##############################
# Outputs directories
##############################
const wmatdir="wmat/"
const gfuncdir="gfunc/"
const axidir = "xifunc/"
const modedir = "xifunc/"

##############################
# Model DF
##############################
const qDF = 6
const σDF = OrbitalElements.σMestelDF(R0,V0,qDF)
const CDF = OrbitalElements.NormConstMestelDF(G,R0,V0,qDF)

const Rin, Rout, Rmax = 1., 11.5, 20.   # Tapering radii
const ξDF = 1.0                         # Self-gravity fraction
const μDF, νDF = 5, 4                   # Tapering exponants

const dfname = "Zang_q_"*string(qDF)*"_xi_"*string(ξDF)*"_mu_"*string(μDF)*"_nu_"*string(νDF)

const DF(E::Float64,L::Float64) = ξDF * OrbitalElements.ZangDF(E,L,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)

const ndFdJ(n1::Int64,n2::Int64,
            E::Float64,L::Float64,
            ndotΩ::Float64)   = ξDF * OrbitalElements.ZangndDFdJ(n1,n2,E,L,ndotΩ,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)

#####
# Parameters
#####
# Radii for frequency truncations
const rmin = 0.1
const rmax = 100.

const Orbitalparams = OrbitalElements.OrbitalParameters(;Ω₀=Ω₀,rmin=rmin,rmax=rmax)


const Ku = 200           # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

const Kv = 201    # number of allocations is directly proportional to this
const Kw = 202    # number of allocations is insensitive to this (also time, largely?

const lharmonic = 2
const n1max = 1  # maximum number of radial resonances to consider


####
const nradial = basis.nmax
const KuTruncation=10000

const VMAPN = 2
const VERBOSE = 1
const OVERWRITE = false

####


const ADAPTIVEKW = true

params = LinearResponse.LinearParameters(;Orbitalparams=Orbitalparams,
                                                Ku=Ku,Kv=Kv,Kw=Kw,
                                                modelname=modelname,dfname=dfname,
                                                wmatdir=wmatdir,gfuncdir=gfuncdir,axidir=axidir,modedir=modedir,
                                                lharmonic=lharmonic,n1max=n1max,nradial=nradial,
                                                KuTruncation=KuTruncation,
                                                VMAPN=VMAPN,
                                                VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                                ndim=basis.dimension,
                                                nmax=basis.nmax,rbasis=basis.rb,ADAPTIVEKW=ADAPTIVEKW)