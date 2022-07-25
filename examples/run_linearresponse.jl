"""

all steps combined into one: could pause and restart if this was too much.


"""



import CallAResponse
using HDF5

inputfile = "ModelParamMestelZang.jl"

CallAResponse.RunWmat(inputfile)

CallAResponse.RunGfunc(inputfile)

# need to organise the omegalist here
#omglist = [0.1-0.2*im]

# might want to consider trying the neutral modes first?

nbEta = 100                               # Number of eta for which the matrix is computed
Etamin = 0.001#*Omega0                   # Minimum eta value for which the ROI matrix is computed
Etamax = 0.030#*Omega0                   # Maximum eta value for which the ROI matrix is computed
deltaEta = (Etamax-Etamin)/(nbEta-1)     # Getting the step distance of the array in eta
tabEta = collect(Etamin:deltaEta:Etamax) # Table of eta for which the response matrix is computed
tabOmg = 0 .- tabEta .* im               # Table of eta for which the response matrix is computed
tabdetXi = zeros(Float64,nbEta)          # Table to store the value of det[I-Xi].

tabdet = CallAResponse.RunM(inputfile,tabOmg)

for iEta=1:nbEta
    println(tabOmg[iEta]," -- ",tabdet[iEta])
end

# for the minimum, go back and compute the mode shape: replace with the minimum adaptively found
EV,EF,EM = CallAResponse.RunShape(inputfile,0 + 0.026im)

println(EV)
println(EF)
println(EM)

ModeR,ModeShape = CallAResponse.tabShapeMode(inputfile,0.01,15.,100,EM)


#println(tabdet)
#println(tabmev)
