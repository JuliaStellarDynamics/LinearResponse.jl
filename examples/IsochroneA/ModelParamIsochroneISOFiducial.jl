"""
an input file which reproduces the Fouvry & Prunet (2022) damped l=1 mode calculation
driven by run_linearresponseIsochrone_damped.jl
"""

using AstroBasis
using DistributionFunctions
using FiniteHilbertTransform
using HDF5
using LinearResponse
using OrbitalElements
using Plots

# choose a basis for computation of the Fourier-transformed basis elements
G         = 1.     # the gravitational constant
rb        = 20.0   # the scale for the basis elements
lmax,nradial = 1,100  # usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

# choose a model potential
modelname = "IsochroneA"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.AnalyticIsochrone()

dfname = "isotropic"
distributionfunction = IsotropicIsochrone(model)


# output directories
wmatdir    = "wmat/"
gfuncdir   = "gfunc/"
modedir    = "xifunc/"

mkpath(wmatdir)
mkpath(gfuncdir)
mkpath(modedir)


Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max     = 10  # maximum number of radial resonances to consider


VERBOSE   = 2
OVERWRITE = false
VMAPN     = 1
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