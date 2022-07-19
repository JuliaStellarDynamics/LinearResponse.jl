

# we will assume you are running from the 'examples' directory
include("../src/Resonances.jl")

lmax  = 1  # maximum harmonic
n1max = 4  # maximum number of radial resonances to consider (usually lmax)

nbResVec = get_nbResVec(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)
println(tabResVec)
println(nbResVec)

nbResVec2 = get_nbResVec_2d(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec2 = maketabResVec_2d(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)
println(tabResVec2)
println(nbResVec2)


# test some loops using parallelised computing
# https://julialang.org/blog/2019/07/multithreading/
# JULIA_NUM_THREADS=4 julia run_resonances.jl

Threads.@threads for i = 1:nbResVec
    n1,n2 = tabResVec[1,i],tabResVec[2,i]
    println("i = $i on thread $(Threads.threadid()) for n1=$n1 and n2=$n2")
end
