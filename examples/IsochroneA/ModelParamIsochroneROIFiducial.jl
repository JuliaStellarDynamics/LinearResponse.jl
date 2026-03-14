"""
an input file which reproduces the Fouvry & Prunet (2022) unstable l=2 mode calculation
driven by runlinearresponseIsochroneunstable.jl
"""

using AstroBasis
using DistributionFunctions
using FiniteHilbertTransform
using HDF5
using LinearResponse
using OrbitalElements

# Basis
G  = 1.
rb = 5.0
lmax,nradial = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# choose a model potential
modelname = "IsochroneA"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.AnalyticIsochrone()

# Model Distribution Function
dfname = "roi1.0"
distributionfunction = OsipkovMerrittIsochroneEL(1.0,model)


# integration parameters

Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max     = 1  # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

mkpath(wmatdir)
mkpath(gfuncdir)
mkpath(modedir)


VERBOSE   = 2
OVERWRITE = false
VMAPN     = 1 # exponent for v mapping (1 is linear)
ADAPTIVEKW= false
KUTRUNCATION=10000

RMIN = 0.0
RMAX = Inf


OEparams = OrbitalElements.OrbitalParameters(rmin=RMIN,rmax=RMAX,
                                             EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ω₀=frequency_scale(model),Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)