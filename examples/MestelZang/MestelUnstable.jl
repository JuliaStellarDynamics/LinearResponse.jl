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


using OrbitalElements
using AstroBasis
using DistributionFunctions
using FiniteHilbertTransform
using LinearResponse
using HDF5



##############################
# Basis
##############################
const G  = 1.

# Clutton-Brock (1972) basis
const basisname = "CluttonBrock"
const rb = 5.
const lmax,nradial = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
const basis = AstroBasis.CB72Basis(lmax=lmax,nradial=nradial,G=G,rb=rb) 

##############################
# Model Potential
##############################
const modelname = "Mestel"

const R0, V0 = 1.,1.#20., 1.
#model = OrbitalElements.MestelPotential(R0=R0,V0=V0)
model = OrbitalElements.TaperedMestel(R0=R0,V0=V0)

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
const Rin, Rout, Rmax = 1., 11.5, 20.   # Tapering radii
const ξDF = 1.0                         # Self-gravity fraction
const μDF, νDF = 2, 2                  # Tapering exponants

# is this meant to be truncated or not truncated?
distributionfunction = TruncatedZangDisc(model,qDF,νDF,Rin,μDF,Rout,Rmax,ξDF,G)

const σDF = σMestelDF(distributionfunction)
const CDF = NormConstMestelDF(distributionfunction)

const dfname = "Zang_q_"*string(qDF)*"_xi_"*string(ξDF)*"_mu_"*string(μDF)*"_nu_"*string(νDF)


"""
const DF(E::Float64,L::Float64) = ξDF * OrbitalElements.ZangDF(E,L,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)

const ndFdJ(n1::Int64,n2::Int64,
            E::Float64,L::Float64,
            ndotΩ::Float64)   = ξDF * OrbitalElements.ZangndDFdJ(n1,n2,E,L,ndotΩ,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)
"""

#####
# Parameters
#####
# OrbitalElements parameters
const EDGE = 0.01
const TOLECC = 0.01
# Radii for frequency truncations
const rmin = 0.2
const rmax = 20.0

const Orbitalparams = OrbitalElements.OrbitalParameters(;rmin=rmin,rmax=rmax,EDGE=EDGE,TOLECC=TOLECC)


const Ku = 200           # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHT(Ku)

const Kv = 200    # number of allocations is directly proportional to this
const Kw = 200    # number of allocations is insensitive to this (also time, largely?

const VMAPN = 2
const KuTruncation=1000

const lharmonic = 2
const n1max = 2  # maximum number of radial resonances to consider


####

const VERBOSE = 1
const OVERWRITE = false

####


const ADAPTIVEKW = false

params = LinearResponse.LinearParameters(basis;Orbitalparams=Orbitalparams,Ω₀=frequency_scale(model),
                                         Ku=Ku,Kv=Kv,Kw=Kw,
                                         VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW,KuTruncation=KuTruncation,
                                         modelname=modelname,dfname=dfname,
                                         wmatdir=wmatdir,gfuncdir=gfuncdir,axidir=axidir,modedir=modedir,
                                         OVERWRITE=OVERWRITE,
                                         lharmonic=lharmonic,n1max=n1max,
                                         VERBOSE=VERBOSE)