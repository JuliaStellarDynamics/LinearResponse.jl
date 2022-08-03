


import CallAResponse
using HDF5

inputfile = "ModelParamIsochrone_damped.jl"

CallAResponse.RunWmat(inputfile)
#CallAResponse.RunWmatIsochrone(inputfile)

CallAResponse.RunGfunc(inputfile)
#CallAResponse.RunGfuncIsochrone(inputfile)

# need to organise the omegalist here
#omglist = [0.1-0.2*im]

nOmega = 52#501 # Number of omega0 for which the response matrix is computed. ATTENTION, we choose this number to have a nice step distance
Omegamin = 0.0#*Omega0 # Minimum omega value for which the response matrix is computed. ATTENTION, we use the symmetry of the response matrix, so that we only search for positive real(omega).
Omegamax = 0.05#1*Omega0 # Maximum omega value for which the response matrix is computed
tabOmega = collect(range(Omegamin,Omegamax,length=nOmega)) # Table of omega for which the response matrix is computed

# Evaluate the dispersion function even in eta=0.0 using the `damped' evaluation,
# but this should not really change the plot that we are making in the lower-half of the complex plane.
nEta = 51#501 # Number of eta for which the response matrix is computed. ATTENTION, we choose this number to have a nice step distance
Etamin = -0.005#-0.06*Omega0#-0.005*Omega0 # Minimum eta value for which the response matrix is computed. ATTENTION, it is a negative number
Etamax = 0.0#-0.0*Omega0 # Maximum eta value for which the response matrix is computed. ATTENTION, it is a negative number
tabEta = collect(range(Etamin,Etamax,length=nEta)) # Table of eta for which the response matrix is computed
#####
nomega = nOmega*nEta # Total number of complex frequencies for which the response matrix is computed
tabomega = zeros(Complex{Float64},nomega) # Table of the complex frequencies for which the response matrix is computed

##################################################
# Function that initialises tabomega
##################################################
function tabomega!()
    icount = 1 # Initialising the counter
    #####
    for iOmega=1:nOmega # Loop over the real part of the frequency
        for iEta=1:nEta # Loop over the complex part of the frequency
            tabomega[icount] = tabOmega[iOmega] + im*tabEta[iEta] # Filling the current value of the complex frequency
            #####
            icount += 1 # Updating the counter
        end
    end
end
##################################################
tabomega!() # Preparing the array of complex frequencies for which the response matrix is computed

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
