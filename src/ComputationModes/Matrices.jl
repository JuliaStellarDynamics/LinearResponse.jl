"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

"""
@TO DESCRIBE
"""
function RunMatrices(ωlist::Array{ComplexF64},
                     FHT::FiniteHilbertTransform.FHTtype,
                     Parameters::ResponseParameters)

    # Preparinng computations of the response matrices
    tabMlist, tabaMcoef, tabωminωmax, FHTlist = PrepareM(Threads.nthreads(),FHT,Parameters)

    # how many omega values are we computing?
    nω = length(ωlist)
    # allocate containers for the matrices
    tabRMreal = zeros(Float64,Parameters.nradial,Parameters.nradial,nω)
    tabRMimag = zeros(Float64,Parameters.nradial,Parameters.nradial,nω)

    (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.RunMatrices: computing $nω frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nω

        k = Threads.threadid()

        if (i==2) && (Parameters.VERBOSE>0) # skip the first in case there is compile time built in
            @time tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],Parameters)
        else
            tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],Parameters)
        end

        for q = 1:Parameters.nradial
            for p = 1:Parameters.nradial
                tabRMreal[p,q,i] = real(tabMlist[k][p,q])
                tabRMimag[p,q,i] = imag(tabMlist[k][p,q])
            end
        end
    end

    h5open(MatFilename(Parameters), "w") do file
        # ω grid
        write(file,"omega",real(ωlist))
        write(file,"eta",imag(ωlist))
        # Matrices
        write(file,"MatricesR",tabRMreal)
        write(file,"MatricesI",tabRMimag)
        # Parameters
        WriteParameters(file,Parameters)
    end

    return tabRMreal, tabRMimag
end