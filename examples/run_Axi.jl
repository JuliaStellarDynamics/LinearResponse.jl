import CallAResponse

inputfile = "ModelParam.jl"

omglist = [0.1-0.2*im]
tabdet, tabmev = CallAResponse.run_M(inputfile,omglist)

println(tabdet)
println(tabmev)