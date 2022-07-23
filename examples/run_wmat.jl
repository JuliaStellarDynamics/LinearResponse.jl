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


import CallAResponse

inputfile = "ModelParam.jl"

CallAResponse.runWmat(inputfile)
