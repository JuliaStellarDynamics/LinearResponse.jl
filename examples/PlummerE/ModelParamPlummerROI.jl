"""
an example input file for running all steps in estimating the Linear Response for a given model

Must include:
-ψ
-dψ
-d2ψ
-basis
-ndim
-nradial
-ndFdJ


"""

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5

#####
# Basis
#####
G  = 1.
rb = 4.0
lmax,nmax = 2,20 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73BasisCreate(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim = basis.dimension
nradial = basis.nmax

#####
# Model Potential
#####


const modelname = "PlummerE"

const bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψPlummer(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

rmin = 1.0e-5
rmax = 1.0e5


dfname = "roi1.0"
dfname = "roi0.75"
#dfname = "roi5.0"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=0.75)

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# integration parameters

Ku = 202    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)


lharmonic = 2
n1max = 2  # maximum number of radial resonances to consider

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 101
Omegamin = -0.1
Omegamax = 0.1
nEta     = 100
Etamin   = -0.1
Etamax   = 0.4



# output directories
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"


VERBOSE   = 2
OVERWRITE = true
EDGE      = 0.01
ELTOLECC  = 0.0005
VMAPN     = 1 # exponent for v mapping (1 is linear)
ADAPTIVEKW= false

OEparams = OrbitalElements.OrbitsParametersCreate(dψ,d2ψ,Ω₀,rmin=rmin,rmax=rmax,
                                                  EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                                  NINT=OrbitalElements.DEFAULT_NINT,FDIFF=OrbitalElements.DEFAULT_TOL,
                                                  da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                                  ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)

# ResponseParametersCreate(OEparams;Ku,Kv,Kw,modelname,dfname,wmatdir,gfuncdir,modedir,lharmonic,n1max,nradial,KuTruncation,VERBOSE,OVERWRITE,ndim,nmax,rbasis,VMAPN,ADAPTIVEKW)
Parameters = LinearResponse.ResponseParametersCreate(OEparams,Ku=Ku,Kv=Kv,Kw=Kw,
                                                    modelname=modelname,dfname=dfname,
                                                    wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,
                                                    lharmonic=lharmonic,n1max=n1max,nradial=nradial,
                                                    KuTruncation=KuTruncation,
                                                    VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                                    ndim=ndim,
                                                    nmax=basis.nmax,rbasis=basis.rb,VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)

# WARNING : / at the end to check !
