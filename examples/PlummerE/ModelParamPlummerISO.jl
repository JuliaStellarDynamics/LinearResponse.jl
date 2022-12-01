"""
an example input file for running all steps in estimating the Linear Response for a given model

Must include:
-potential
-dpotential
-ddpotential
-basis
-ndim
-nradial
-ndFdJ


"""


import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import CallAResponse
using HDF5

# basis parameters
G  = 1.
rb = 20.0
lmax,nmax = 1,100 # number of basis functions

# automatically set up the basis parameters
basis   = AstroBasis.CB73BasisCreate(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim    = basis.dimension
nradial = basis.nmax

rmin = 1.0e-5
rmax = 1.0e5


# model Potential
modelname = "PlummerE"
bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψPlummer(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

#dfname = "isotropic"

dfname = "roiinf"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=10000.0)

    #return OrbitalElements.plummer_ISO_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG)
    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# integration parameters

Ku = 202    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000

# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)


lharmonic = 1
n1max = 5  # maximum number of radial resonances to consider

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


VERBOSE   = 2
OVERWRITE = true
EDGE      = 0.01
ELTOLECC  = 0.0005
VMAPN     = 1 # exponent for v mapping (1 is linear)


Parameters = CallAResponse.ResponseParametersCreate(dψ,d2ψ,Ku=Ku,Kv=Kv,Kw=Kw,
                                                    modelname=modelname,dfname=dfname,
                                                    wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,
                                                    lharmonic=lharmonic,n1max=n1max,nradial=nradial,
                                                    KuTruncation=KuTruncation,
                                                    VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                                    Ω₀=Ω₀,rmin=rmin,rmax=rmax,
                                                    EDGE=EDGE,ELTOLECC=ELTOLECC,ndim=ndim,
                                                    nmax=basis.nmax,rbasis=basis.rb,VMAPN=VMAPN)




# WARNING : / at the end to check !
