"""
the input file for computing the Fiducial isotropic calculation in the isochrone case
"""


using AstroBasis
using DistributionFunctions
using FiniteHilbertTransform
using HDF5
using LinearResponse
using OrbitalElements
using Plots


# Basis
G  = 1.
rb = 5.0
lmax,nradial = 1,20#100

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

# Model Potential
const modelname = "IsochroneE2"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.NumericalIsochrone()

# rmin = 0.0
# rmax = Inf


dfname = "isotropic"
distributionfunction = IsotropicIsochrone(model)



# integration parameters
Ku = 20
Kv = 20
Kw = 20

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = lmax
n1max = 4  # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

mkpath(wmatdir)
mkpath(gfuncdir)
mkpath(modedir)


VERBOSE   = 1
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

# WARNING : / at the end to check !
