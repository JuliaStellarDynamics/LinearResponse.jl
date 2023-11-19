"""
an example input file for running all steps in estimating the Linear Response for the isotropic Plummer model

"""

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5

# Basis
G  = 1.
rb = 5.0
lmax,nradial = 1,100 # Usually lmax corresponds to the considered harmonics lharmonic

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# Model Potential
const modelname = "PlummerE"
const bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

# Model DF
dfname = "roiinf"
function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=100000.0)

    return OrbitalElements.plummer_ISO_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG)

end



# integration parameters
Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = lmax
n1max = 10  # maximum number of radial resonances to consider

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 40
Omegamin = 0.0
Omegamax = 0.1
nEta     = 40
Etamin   = -0.1
Etamax   = 0.4

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

# WARNING : / at the end to check !
