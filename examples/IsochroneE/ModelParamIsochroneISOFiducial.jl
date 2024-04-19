"""
the input file for computing the Fiducial isotropic calculation in the isochrone case
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
lmax,nradial = 1,100

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

# Model Potential
const modelname = "IsochroneE2"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.NumericalIsochrone()

rmin = 0.0
rmax = Inf


dfname = "isotropic"
distributionfunction = IsotropicIsochrone(model)



# integration parameters
Ku = 20
Kv = 20
Kw = 20

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = 1
n1max = 4  # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 40
Omegamin = 0.0
Omegamax = 0.1
nEta     = 40
Etamin   = -0.1
Etamax   = 0.4

VERBOSE   = 1
OVERWRITE = false
VMAPN     = 1
ADAPTIVEKW= false
KUTRUNCATION=10000

# use almost entirely default parameters
OEparams = OrbitalElements.OrbitalParameters(EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ω₀=OrbitalElements.frequency_scale(model),Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)


# WARNING : / at the end to check !
