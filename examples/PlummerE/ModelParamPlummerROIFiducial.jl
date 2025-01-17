"""
the input file to compute the radial orbit instability in the Plummer sphere, using Fiducial parameters
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
rb = 3.0
lmax,nradial = 1,25 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# Model Potential
const modelname = "PlummerE"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.NumericalPlummer()

# Model Distribution Function
dfname = "roi0.75"
distributionfunction = OsipkovMerrittPlummerEL(0.75,model)



# Linear Response integration parameters
Ku = 150    # number of Legendre integration sample points
Kv = 150    # number of allocations is directly proportional to this
Kw = 150    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# Define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = lmax
n1max = 10  # maximum number of radial resonances to consider
n1max = 1  # the Fiducial value is 10, but in the interest of a quick calculation, we limit ourselves to 1.

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"


VERBOSE   = 2
OVERWRITE = false
VMAPN     = 1
ADAPTIVEKW= false

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
