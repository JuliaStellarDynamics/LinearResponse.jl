
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneROIhires.jl"
include(inputfile)

import LinearResponse
using HDF5

modeRmin = 0.01
modeRmax = 24.0
nmode = 800

# compute the Fourier-transformed basis elements
LinearResponse.RunWmat(ψ,dψ,d2ψ,FHT,basis,Parameters)
println("Done with W matrices.")

function ndFdJplaceholder(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.6)

    Q = OrbitalElements.isochroneQROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochroneSahadDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochronedQdEROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochronedQdLROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL

    F = OrbitalElements.isochroneSahaDF(E,L,Ra,bc,M,astronomicalG)

    smoothingprefactor = 0.01 # zz

    smoothingprefactor = 0.001 # yy
    smoothingprefactor = 0.02 # xx
    smoothingprefactor = 0.005 # ww

    smoothing = exp(-smoothingprefactor/Q)
    dsmoothing = smoothingprefactor * exp(-smoothingprefactor/Q) / (Q^2)

    smoothdFdQ = smoothing*dFdQ + F*dsmoothing

    #res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ
    res = smoothdFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end

minmode = 0.88
#minmode = 1.0
#minmode = 1.4

startinggamma = 0.02
startinggammas = [1.e-5,0.001,0.002,0.005,0.012,0.02,0.03]

# open the file
#open("mode"*string(startinggamma)*"a.txt","w") do io
open("mode_ww_nrad"*string(basis.nradial)*".txt","w") do io

    for rastep=1:60#120

        raval = round((rastep-1)*0.02 + minmode,digits=3)

        dfname = "roi$raval"

        Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ku=Ku,Kv=Kv,Kw=Kw,
                                                     modelname=modelname,dfname=dfname,
                                                     wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                                     lharmonic=lharmonic,n1max=n1max,
                                                     KuTruncation=KuTruncation,
                                                     VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                                     VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)

        # reset the distribution function call
        ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64) = ndFdJplaceholder(n1,n2,E,L,ndotOmega,Ra=raval)

        # compute the G(u) functions
        LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

        # compute the matrix response at each location in tabomega
        tabdet = LinearResponse.RunAXi(FHT,Parameters)


        #startingomg = 0.0 - 0.0000011773114542602149im
        # find a pole by using gradient descent

        for g=1:length(startinggammas)
            startingomg = 0.0 + startinggammas[g]*1im

            POLETOL,DSIZE = 1.e-16,1.e-5
            bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters,POLETOL,DSIZE)

            # print to file: raval, real omega, imag omega, real determinant, imag determinant, launching gamma
            println(io,raval," ",real(bestomg)," ",imag(bestomg)," ",real(detval)," ",imag(detval)," ",startinggammas[g])

            #bestomg = 0.0 + 0.02271406012170436im
            println("The $raval zero-crossing frequency is $bestomg, with determinant $detval.")

            if ((isfinite(bestomg)) & (abs(real(detval)) < 1.0))
                # for the minimum, go back and compute the mode shape
                EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

                ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters,output=false)

                h5open(Parameters.modedir*"ModeShape_"*string(startinggammas[g])*"_"*LinearResponse.OutputsEnd(Parameters), "w") do file
                    write(file,"ModeRadius",ModeRadius)         # write tabRMode to file
                    write(file,"ModePotentialShape",ModePotentialShape) # write tabShapeMode to file
                    write(file,"ModeDensityShape",ModeDensityShape) # write tabShapeMode to file
                    # Parameters
                    LinearResponse.WriteParameters(file,Parameters)
                end
            end
        end

    end # end loop over ra vals

end # file closure
