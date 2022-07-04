"""
options to run:

julia run_wmat.jl
# 212.325749 seconds (2.34 G allocations: 50.334 GiB, 3.81% gc time, 2.25% compilation time)

julia --compile=min run_wmat.jl
# This one is essentially impossible: multiple hours!

julia --optimize=3 run_wmat.jl
# identical to default

julia --optimize=0 run_wmat.jl
# 312.212810 seconds (2.55 G allocations: 53.461 GiB, 2.80% gc time, 1.19% compilation time)

julia --optimize=3 --inline=no run_wmat.jl
1136.466458 seconds (6.94 G allocations: 141.343 GiB, 1.44% gc time, 0.73% compilation time)

julia --compile=all run_wmat.jl
# 423.578417 seconds (2.34 G allocations: 50.334 GiB, 3.91% gc time, 2.05% compilation time)

@IMPROVE, where should we parallelize this?

"""


import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


# we will assume you are running from the 'examples' directory
include("../src/WMat.jl")
basedir=""

# set up the AstroBasis call
rb,G = 5.,1.
lmax,nmax = 2,10
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
AstroBasis.fill_prefactors!(basis)

# construct the table of needed resonance vectors


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
