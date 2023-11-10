
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamIsochroneROI.jl"
include(inputfile)

import LinearResponse
using HDF5

modeRmin = 0.01
modeRmax = 18.0
nmode = 500

# compute the Fourier-transformed basis elements
#LinearResponse.RunWmatIsochrone(FHT,bc,M,G,basis,Parameters)

function ndFdJplaceholder(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.6)

    Q = OrbitalElements.isochroneQROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochroneSahadDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochronedQdEROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochronedQdLROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end

minmode = 0.87
#minmode = 1.0
#minmode = 1.4

# open the file
open("mode0p003.txt","w") do io

    for rastep=1:40

        # find a pole by using gradient descent
        startingomg = 0.0 + 0.003im

        raval = round((rastep-1)*0.01 + minmode,digits=3)

        if raval==1.34
            continue
        end

        dfname = "roi$raval"

        Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ku=Ku,Kv=Kv,Kw=Kw,
                                                     modelname=modelname,dfname=dfname,
                                                     wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                                     lharmonic=lharmonic,n1max=n1max,
                                                     KuTruncation=KuTruncation,
                                                     VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                                     VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)

        #println(Parameters.dfname)

        if rastep==35
            continue
        end

        if rastep==65
            continue
        end

        # reset the distribution function call
        ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64) = ndFdJplaceholder(n1,n2,E,L,ndotOmega,Ra=raval)

        # compute the G(u) functions
        LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

        # compute the matrix response at each location in tabomega
        tabdet = LinearResponse.RunAXi(FHT,Parameters)


        #startingomg = 0.0 - 0.0000011773114542602149im
        bestomg = LinearResponse.FindPole(startingomg,FHT,Parameters,1.e-12)

        println(io,raval," ",real(bestomg)," ",imag(bestomg))

        #bestomg = 0.0 + 0.02271406012170436im
        println("The $raval zero-crossing frequency is $bestomg.")

        if isfinite(bestomg)
            # for the minimum, go back and compute the mode shape
            EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

            ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)
        end

    end # end loop over ra vals

end # file closure
