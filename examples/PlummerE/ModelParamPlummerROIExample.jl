"""
the input file to compute the radial orbit instability in the Plummer sphere, using Fiducial parameters
"""

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5

# Basis
G  = 1.
rb = 3.0
lmax,nradial = 2,5 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# Model Potential
const modelname = "PlummerE"
const bc, M = 1.,1. # G is defined above: must agree with basis!
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

# Model Distribution Function
dfname = "roi0.75"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=0.75)

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# Linear Response integration parameters
Ku = 12    # number of Legendre integration sample points
Kv = 20    # number of allocations is directly proportional to this
Kw = 20    # number of allocations is insensitive to this (also time, largely)?


# Define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max = 1  # the Fiducial value is 10, but in the interest of a quick calculation, we limit ourselves to 1.

# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 120
Omegamin = -0.1
Omegamax = 0.1
nEta     = 120
Etamin   = -0.05
Etamax   = 0.1


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
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)

# WARNING : / at the end to check !
