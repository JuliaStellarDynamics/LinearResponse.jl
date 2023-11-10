
# to mimic Fouvry & Prunet exactly
inputfile = "ModelParamPlummerROIhires.jl"
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

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


println(ndFdJplaceholder(-1,2,-1.,1.,1.))

minmode = 0.75

startinggammas = [1.e-5,0.01,0.02,0.05]

# open the file
#open("mode"*string(startinggamma)*"a.txt","w") do io
open("mode_adjust2b_nrad"*string(basis.nradial)*".txt","w") do io

    for rastep=1:60

        raval = round((rastep-1)*0.005 + minmode,digits=4)

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
            bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)#,POLETOL,DSIZE)

            # print to file: raval, real omega, imag omega, real determinant, imag determinant, launching gamma
            println(io,raval," ",real(bestomg)," ",imag(bestomg)," ",real(detval)," ",imag(detval)," ",startinggammas[g])

            #bestomg = 0.0 + 0.02271406012170436im
            println("The $raval zero-crossing frequency is $bestomg, with determinant $detval.")

            if ((isfinite(bestomg)) & (abs(real(detval)) < 1.0))
                # for the minimum, go back and compute the mode shape
                EV,EM = LinearResponse.ComputeModeTables(bestomg,FHT,Parameters)

                ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters,output=false)
                println(Parameters.modedir*"ModeShape_"*string(startinggammas[g])*"_"*LinearResponse.OutputsEnd(Parameters))
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
