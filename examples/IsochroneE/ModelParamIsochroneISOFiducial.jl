"""
the input file for computing the Fiducial isotropic calculation in the isochrone case
"""


import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5


# Basis
G  = 1.
rb = 5.0
lmax,nradial = 1,100

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

# Model Potential
modelname = "IsochroneE"
bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)


# define the distribution function and give it a name
dfname = "isotropic"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    dFdE = OrbitalElements.isochroneisotropicdDFdE(E,bc,M,astronomicalG)
    res = ndotOmega*dFdE
    return res
end


# integration parameters
Ku = 200
Kv = 200
Kw = 200

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = 1
n1max = 10  # maximum number of radial resonances to consider

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
OEparams = OrbitalElements.OrbitalParameters(Ω₀=Ω₀,
                                             EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             KuTruncation=KUTRUNCATION,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)


# WARNING : / at the end to check !
