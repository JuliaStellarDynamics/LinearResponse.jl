import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5
include("/Users/mpetersen/CodeHold/JuliaCallAResponse/src/WMat.jl")

basedir="/Users/mpetersen/CodeHold/JuliaCallAResponse/examples/"

# set up the AstroBasis call
rb,G = 5.,1.
Ulnp,DF = AstroBasis.read_and_fill_prefactors(2,10,rb,G)

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
n1 = -2
n2 = -2
NstepsWMat = 50

for n1 in [-2,-1,0,1,2]
    for n2 in [-2,0,2]
        if (n1==0) & (n2==0)
            continue
        end
        println(n1," ",n2)
        tabWMat,tabaMat,tabeMat = make_wmat(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,nradial,lharmonic,Ulnp,Omega0,rb,NstepsWMat)
        # now save
        h5open(basedir*"wmat/wmat_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*".h5", "w") do file
            write(file, "wmat",tabWMat)
            write(file, "amat",tabaMat)
            write(file, "emat",tabeMat)
        end
    end
end
