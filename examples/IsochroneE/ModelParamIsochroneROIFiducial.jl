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
rb = 5.0
lmax,nradial = 2,100 # Usually lmax corresponds to the considered harmonics lharmonic

# CB73Basis([name, dimension, lmax, nradial, G, rb, filename])
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)


# Model Potential
modelname = "IsochroneE"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)


rmin = 0.0
rmax = Inf


dfname = "roi1.0"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)

    Q = OrbitalElements.isochroneQROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochroneSahadDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochronedQdEROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochronedQdLROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end


# integration parameters

Ku = 200    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000


# define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


lharmonic = 2
n1max = 1  # maximum number of radial resonances to consider

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 11
Omegamin = -0.02
Omegamax = 0.02
nEta     = 10
Etamin   = 0.001
Etamax   = 0.04


# output directories
wmatdir  = "wmat/"
#gfuncdir = "/Volumes/External1/P23RegressionTests/"#gfunc/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"


VERBOSE   = 2
OVERWRITE = true
EDGE      = 0.01
ELTOLECC  = 0.01
VMAPN     = 1
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


# WARNING : / at the end to check !
