


import CallAResponse
using HDF5

inputfile = "ModelParamIsochrone_damped.jl"

#CallAResponse.RunWmat(inputfile)
#CallAResponse.RunWmatIsochrone(inputfile)

#CallAResponse.RunGfunc(inputfile)
#CallAResponse.RunGfuncIsochrone(inputfile)

# need to organise the omegalist here

# need this to get the parameters... is this the best fix?
include(inputfile)

tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

#=
tabdet = CallAResponse.RunM(inputfile,tabomega)

# make square versions of everything
omegatable = zeros(nOmega,nEta)
etatable = zeros(nOmega,nEta)
dettable = zeros(nOmega,nEta)

for iOmega=1:nOmega # Loop over the real part of the frequency
    for iEta=1:nEta
        combo = (iOmega-1)*(nOmega-1) + iEta

        omegatable[iOmega,iEta] = real(tabomega[combo])#tabOmega[iOmega]
        etatable[iOmega,iEta] = imag(tabomega[combo])#tabEta[iEta]
        #println("iO=$iOmega,iE=$iEta,combo=$combo,O=$(tabOmega[iOmega]),E=$(tabEta[iEta])")
        #println("iO=$iOmega,iE=$iEta,combo=$combo,O=$(real(tabomega[combo])),E=$(imag(tabomega[combo]))")
        dettable[iOmega,iEta] = tabdet[combo]
    end
end

# now need to print these
#h5open(det_filename(modedir,modelname,lharmonic,n1max,K_u), "w") do file
detname = "xifunc/Determinant_test.h5"

h5open(detname, "w") do file
    write(file,"omega",real(tabomega))
    write(file,"eta",imag(tabomega))
    write(file,"det",tabdet)
end

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
=#
