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
using ArgParse # To parse command-line arguments

#####
# Input file parsing and check
#####
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--inputfile"
    help = "File with test case specifications."
    arg_type = String
    default = "./test.jl"
end
parsed_args = parse_args(tabargs)

input = parsed_args["inputfile"] 
include(input)

if !( (@isdefined G) && (@isdefined rb) && (@isdefined basis) 
    && (@isdefined modelname) 
    && (@isdefined potential) && (@isdefined dpotential) && (@isdefined ddpotential) 
    && (@isdefined K_u) && (@isdefined K_v) && (@isdefined NstepsWMat) 
    && (@isdefined lharmonic) && (@isdefined n1max) 
    && (@isdefined wmatdir) )
    
    error("Definitions missing among G, rb, basis, 
            modelname, potential, dpotential, ddpotential, 
            K_u, K_v, NstepsWMat, lharmonic, n1max, wmatdir")
end
if last(wmatdir) != '/'
    error(" '/' should be included at the end of wmatdir")
end

#####
# Computation
#####

# we will assume you are running from the 'examples' directory
include("../src/WMat.jl")          # fully adaptive WMat calculations (production)
#include("../src/WMatIsochrone.jl") # for isochrone-specific WMat calculations (testing purposes)
include("../src/Resonances.jl")    # for resonances helpers


#####
# Bases prep.
#####
AstroBasis.fill_prefactors!(basis)
bases=[deepcopy(basis) for k=1:Threads.nthreads()]



#####
# Legendre integration prep.
#####
tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

#####
# Construct the table of needed resonance vectors
#####
nbResVec = get_nbResVec_2d(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec_2d(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)

println(nbResVec)

Threads.@threads for i = 1:nbResVec
    k = Threads.threadid()
    n1,n2 = tabResVec[1,i],tabResVec[2,i]
    println(n1," ",n2)
    @time tabWMat,tabaMat,tabeMat = make_wmat(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,lharmonic,bases[k],Omega0)
    #make_wmat_isochrone(potential,dpotential,ddpotential,n1,n2,tabuGLquad,K_v,lharmonic,basis,Omega0)
    #println(sum(tabWMat))
    # now save
    h5open(wmatdir*"wmat_"*string(modelname)*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rb)*".h5", "w") do file
        write(file, "wmat",tabWMat)
        write(file, "amat",tabaMat)
        write(file, "emat",tabeMat)
    end
end