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


#####
# Basis
#####
G         = 1.
rb        = 15.0
lmax,nmax = 1,30
basis     = AstroBasis.CB73BasisCreate(lmax=lmax, nmax=nmax,G=G,rb=rb)
ndim      = basis.dimension
nradial   = basis.nmax

#########
# Model Potential
#########
modelname = "IsochroneE"

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)
rmin = 1.0e-5
rmax = 1.0e5


#########
# Distribution Function
#########
dfname = "isotropic"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG) # Current value of dF/E. ATTENTION, the DF is assumed to be isotropic
    res = ndotOmega*dFdE # Current value of n.dF/dJ. ATTENTION, the DF is assumed to be isotropic
    #####
    return res
end


#########
# Integration parameters
#########
Ku = 203    # number of Legendre integration sample points
Kv = 200    # number of allocations is directly proportional to this
Kw = 200    # number of allocations is insensitive to this (also time, largely)?
KuTruncation = 10000 # if limiting Ku for sum, specify here

FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

#########
# Considered resonance parameters
#########
lharmonic = 1
n1max     = 2  # maximum number of radial resonances to consider

#########
# Frequencies to probe
#########
nOmega   = 51
Omegamin = 0.0
Omegamax = 0.05
nEta     = 50
Etamin   = -0.005
Etamax   = 0.0


#########
# output directories
#########
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

#########
# control parameters
#########
VERBOSE   = 2
OVERWRITE = false
EDGE      = 0.01
ELTOLECC  = 0.0005
VMAPN     = 2 # exponent for v mapping (1 is linear)
ADAPTIVEKW= true

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
