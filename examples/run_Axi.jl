import OrbitalElements
import PerturbPlasma
using HDF5
using LinearAlgebra

include("/Users/mpetersen/CodeHold/JuliaCallAResponse/src/Xi.jl")
include("/Users/mpetersen/CodeHold/JuliaCallAResponse/src/Basis.jl")

basedir="/Users/mpetersen/CodeHold/JuliaCallAResponse/examples/"

bc, M, G = 1.,1. ,1.
potential   = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential  = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)

# bring in Legendre integration prefactors
K_u = 200
# get all weights
tabuGLquadtmp,tabwGLquadtmp,tabINVcGLquadtmp,tabPGLquadtmp = PerturbPlasma.tabGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)
tabwGLquad = reshape(tabwGLquadtmp,K_u,1)
tabINVcGLquad = reshape(tabINVcGLquadtmp,K_u,1)
tabPGLquad = reshape(tabPGLquadtmp,K_u,1)

lmax=2
n1max=2
nbResVec = get_nbResVec(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)

# make the (np,nq) vectors that we need to evaluate
tab_npnq = makeTabnpnq(nradial)

# initialise the array for memory saving
tabaXi = [[zeros(Float64,K_u) for np=1:nradial,nq=1:nradial] for ResVec=1:nbResVec]

# make the Xi coefficients
makeXiCoefficients!(tabaXi,tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,basedir)

# once we have tabaXi, just go for it!
PARALLEL=true
struct_tabLeg = PerturbPlasma.initialize_struct_tabLeg(K_u,PARALLEL)

# test an example
#w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)
#omg = 0.0 + 0.024im
#omg_nodim = omg/Omega0 # Dimensionless frequency rescaled by Omega0
#n1,n2 = -1,2
#varpi = OrbitalElements.get_varpi(omg_nodim,n1,n2,dpotential,ddpotential)
#PerturbPlasma.get_tabLeg!(varpi,K_u,struct_tabLeg[1],"Unstable")
#struct_tabLeg[1].tabDLeg

# now prepare to do the integration by populating for whichever mode...
tabXi = zeros(Complex{Float64},nradial,nradial)
omg = 0.0 + 0.024im
tabXi!(omg,tabXi,tabaXi,tabResVec,tab_npnq,struct_tabLeg[1],dpotential,ddpotential,"Unstable")

# wrapper bringing it all together
detXi(omg,tabXi,tabaXi,tabResVec,tab_npnq,struct_tabLeg[1],dpotential,ddpotential,"Unstable",Omega0)
