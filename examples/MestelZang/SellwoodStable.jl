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
import LinearResponse

using HDF5


##############################
# Basis
##############################
const G  = 1.

# Clutton-Brock (1972) basis
const basisname = "CluttonBrock"
const rb = 5.
const lmax,nradial = 2,10 # Usually lmax corresponds to the considered harmonics lharmonic
const basis = AstroBasis.CB72Basis(lmax=lmax,nradial=nradial,G=G,rb=rb) 

# # Kalnajs (1976) basis
# basisname = "Kalnajs"
# rb, kKA = 5., 7
# lmax,nmax = 2,7
# basis = AstroBasis.K76Basis(lmax=lmax,nradial=nradial,G=G,rb=rb,kKA=kKA)

##############################
# Model Potential
##############################
const modelname = "Mestelv1"

const R0, V0 = 20., 1.
#const ψ(r::Float64)   = OrbitalElements.ψMestel(r,R0,V0)
#const dψ(r::Float64)  = OrbitalElements.dψMestel(r,R0,V0)
#const d2ψ(r::Float64) = OrbitalElements.d2ψMestel(r,R0,V0)
#const Ω₀ = OrbitalElements.Ω₀Mestel(R0,V0)
model = OrbitalElements.MestelPotential(R0,V0)
println(OrbitalElements.Ω₀(model))


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
const qDF = 11.4
const σDF = OrbitalElements.σMestelDF(R0,V0,qDF)
const CDF = OrbitalElements.NormConstMestelDF(G,R0,V0,qDF)

const Rin, Rout, Rmax = 1., 11.5, 20.   # Tapering radii
const ξDF = 0.5                         # Self-gravity fraction
const μDF, νDF = 5, 4                   # Tapering exponants

const dfname = "Zang_q_"*string(qDF)*"_xi_"*string(ξDF)*"_mu_"*string(μDF)*"_nu_"*string(νDF)

const ndFdJ(n1::Int64,n2::Int64,
            E::Float64,L::Float64,
            ndotΩ::Float64) = ξDF * OrbitalElements.ZangndDFdJ(n1,n2,E,L,ndotΩ,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)


#####
# Parameters
#####
# OrbitalElements parameters
const EDGE = 0.01
const TOLECC = 0.01
# Radii for frequency truncations
const rmin = 0.1
const rmax = 100.0

const Orbitalparams = OrbitalElements.OrbitalParameters(;Ω₀=OrbitalElements.Ω₀(model),rmin=rmin,rmax=rmax,EDGE=EDGE,TOLECC=TOLECC)


const Ku = 100           # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHT(Ku)

const Kv = 101    # number of allocations is directly proportional to this
const Kw = 102    # number of allocations is insensitive to this (also time, largely?

const VMAPN = 2
const KuTruncation=1000

const lharmonic = 2
const n1max = 1  # maximum number of radial resonances to consider


####

const VERBOSE = 1
const OVERWRITE = true

####


const ADAPTIVEKW = true

params = LinearResponse.LinearParameters(basis;Orbitalparams=Orbitalparams,
                                         Ku=Ku,Kv=Kv,Kw=Kw,
                                         VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW,KuTruncation=KuTruncation,
                                         modelname=modelname,dfname=dfname,
                                         wmatdir=wmatdir,gfuncdir=gfuncdir,axidir=axidir,modedir=modedir,
                                         OVERWRITE=OVERWRITE,
                                         lharmonic=lharmonic,n1max=n1max,
                                         VERBOSE=VERBOSE)