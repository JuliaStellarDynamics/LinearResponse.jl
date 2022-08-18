


import CallAResponse
using HDF5

inputfile = "ModelParamIsochrone_damped.jl"

# compute the Fourier-transformed basis elements
CallAResponse.RunWmat(inputfile)

# compute the G(u) functions
CallAResponse.RunGfunc(inputfile)

# need to organise the omegalist here


# now need to print these
#h5open(det_filename(modedir,modelname,lharmonic,n1max,K_u), "w") do file
detname = "xifunc/Determinant_test.h5"

#=


=#

#=
h5open(detname, "w") do file
    write(file,"omega",omegatable)
    write(file,"eta",etatable)
    write(file,"det",dettable)
end
=#

#=
for iEta=1:nbEta
    println(tabOmg[iEta]," -- ",tabdet[iEta])
end

# for the minimum, go back and compute the mode shape
EV,EF,EM = CallAResponse.RunShape(inputfile,0 + 0.026im)

println(EV)
println(EF)
println(EM)

ModeR,ModeShape = CallAResponse.tabShapeMode(inputfile,0.01,15.,100,EM)
=#
