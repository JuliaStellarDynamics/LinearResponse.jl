import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5

# options to run:
# julia --compile=min run_wmat.jl
# julia --optimize=0 run_wmat.jl
# julia --compile=min --inline=no run_wmat.jl

# we will assume you are running from the 'examples' directory
include("../src/WMat.jl")
basedir=""

# set up the AstroBasis call
rb,G = 5.,1.
lmax,nmax = 2,10
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
AstroBasis.fill_prefactors!(basis)

# bring in Legendre integration prefactors
K_u = 200
tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

lharmonic = 2

bc, M, G = 1.,1. ,1.
potential   = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential  = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)

K_v = 50
nradial=5
NstepsWMat = 50

#for n1 in [-2,-1,0,1,2]
#    for n2 in [-2,0,2]
for n1 in [0]
    for n2 in [2]
        if (n1==0) & (n2==0)
            continue
        end
        println(n1," ",n2)
        @time tabWMat,tabaMat,tabeMat = make_wmat(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,lharmonic,basis,Omega0)
        # now save
        h5open(basedir*"wmat/wmat_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*".h5", "w") do file
            write(file, "wmat",tabWMat)
            write(file, "amat",tabaMat)
            write(file, "emat",tabeMat)
        end
    end
end
