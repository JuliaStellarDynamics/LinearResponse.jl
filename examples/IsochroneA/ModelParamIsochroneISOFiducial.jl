"""
an input file which reproduces the Fouvry & Prunet (2022) damped l=1 mode calculation
driven by run_linearresponseIsochrone_damped.jl
"""

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5

# choose a basis for computation of the Fourier-transformed basis elements
G         = 1.     # the gravitational constant
rb        = 20.0   # the scale for the basis elements
lmax,nradial = 1,100  # usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

# choose a model potential
modelname = "IsochroneA"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)



# choose a distribution function for G(u) calculation
dfname = "isotropic"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    dFdE = OrbitalElements.isochroneisotropicdDFdE(E,bc,M,astronomicalG)
    res = ndotOmega*dFdE
    return res
end


# output directories
wmatdir    = "wmat/"
gfuncdir   = "gfunc/"
modedir    = "xifunc/"


Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = 1
n1max     = 10  # maximum number of radial resonances to consider

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

VERBOSE   = 2
OVERWRITE = false
VMAPN     = 1
ADAPTIVEKW= false

OEparams = OrbitalElements.OrbitalParameters(Ω₀=Ω₀,
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
