import CallAResponse

inputfile = "ModelParam.jl"

# need to organise the omegalist here
omglist = [0.1-0.2*im]
tabdet, tabmev = CallAResponse.runM(inputfile,omglist)

println(tabdet)
println(tabmev)
