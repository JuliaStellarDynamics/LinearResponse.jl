


import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5


#####
# Basis
#####
G  = 1.
rb = 20.0
lmax,nradial = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
@time basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# Model Potential
modelname = "IsochroneEt"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)


rmin = 0.0
rmax = Inf
dfname = "roi1.0"


# integration parameters

Ku = 205    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000


# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = 2
n1max = 0  # maximum number of radial resonances to consider

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 51
Omegamin = -0.02
Omegamax = 0.02
nEta     = 50
Etamin   = 0.001
Etamax   = 0.04


# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"


VERBOSE   = 1
OVERWRITE = false
EDGE      = 0.01
ELTOLECC  = 0.0005
VMAPN     = 2 # exponent for v mapping (1 is linear)
ADAPTIVEKW= false

OEparams = OrbitalElements.OrbitalParameters(Ω₀=Ω₀,rmin=rmin,rmax=rmax,
                                             EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             KuTruncation=KuTruncation,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)



LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)
