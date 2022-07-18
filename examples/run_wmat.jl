"""
options to run:

julia run_wmat.jl
# 212.325749 seconds (2.34 G allocations: 50.334 GiB, 3.81% gc time, 2.25% compilation time)
# after converting definitions to const, running on the new machine:
# 73.179156 seconds (2.34 G allocations: 50.337 GiB, 6.72% gc time, 2.05% compilation time)

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

julia run_wmat.jl #with isochrone specific
2.012504 seconds (58.03 M allocations: 1.039 GiB, 5.09% gc time, 35.71% compilation time)


JULIA_NUM_THREADS=4 julia run_wmat.jl

@IMPROVE, where should we parallelize this? (while looping through the resonances)

"""


import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


# we will assume you are running from the 'examples' directory
include("../src/WMat.jl")          # fully adaptive WMat calculations (production)
#include("../src/WMatIsochrone.jl") # for isochrone-specific WMat calculations (testing purposes)
include("../src/Resonances.jl")    # for resonances helpers
basedir=""

# set up the AstroBasis call
const G  = 1.
const rb = 10.
lmax,nmax = 2,10
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
AstroBasis.fill_prefactors!(basis)

# construct the table of needed resonance vectors


# bring in Legendre integration prefactors and sample points
const K_u = 150
tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

const lharmonic = 2

#=
#const modelname = "isochrone"
const modelname = "isochroneE"

const bc, M = 1.,1.
potential(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
=#

const modelname = "PlummerE"

const bc, M = 1.,1.
potential(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.plummer_Omega0(bc,M,G)

const K_v        = 100  # number of allocations is directly proportional to this
const NstepsWMat = 100   # number of allocations is insensitive to this (also time, largely?)
# const nradial?

lmax  = 2  # maximum harmonic
n1max = 4  # maximum number of radial resonances to consider

nbResVec = get_nbResVec(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)

println(nbResVec)

Threads.@threads for i = 1:nbResVec
    n1,n2 = tabResVec[1,i],tabResVec[2,i]
    println(n1," ",n2)
    @time tabWMat,tabaMat,tabeMat = make_wmat(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,lharmonic,basis,Omega0)
    #make_wmat_isochrone(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,lharmonic,basis,Omega0)
    #println(sum(tabWMat))
    # now save
    h5open(basedir*"wmat/wmat_"*string(modelname)*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rb)*".h5", "w") do file
        write(file, "wmat",tabWMat)
        write(file, "amat",tabaMat)
        write(file, "emat",tabeMat)
    end
end
