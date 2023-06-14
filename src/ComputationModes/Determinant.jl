"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""


"""
    detXi(IMat,tabM[,ξ])

determinant of the susceptibility matrix I - ξ*M(ω) 
for known M(ω) and active fraction ξ (default 1).
"""
function detXi(IMat::AbstractMatrix{ComplexF64},
               tabM::AbstractMatrix{ComplexF64};
               ξ::Float64=1.0)::ComplexF64

    # Computing the determinant of (I-M).
    # ATTENTION, we tell julia that the matrix is symmetric
    val = det(Symmetric(IMat-ξ*tabM))

    # only save the real portion
    return val # Output
end


"""
@TO DESCRIBE
"""
function RunDeterminant(ωlist::Array{ComplexF64},
                        FHT::FiniteHilbertTransform.AbstractFHT,
                        params::LinearParameters;
                        ξ::Float64=1.0)

    # Preparinng computations of the response matrices
    tabMlist, tabaMcoef, tabωminωmax, FHTlist = PrepareM(Threads.nthreads(),FHT,params)
    IMatlist = makeIMat(params.nradial,Threads.nthreads())

    # how many omega values are we computing?
    nω = length(ωlist)
    # allocate containers for determinant and min eigenvalue
    tabdetXi = zeros(ComplexF64,nω)

    (params.VERBOSE > 0) && println("LinearResponse.Xi.RunDeterminant: computing $nω frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nω

        k = Threads.threadid()

        if (i==2) && (params.VERBOSE>0) # skip the first in case there is compile time built in
            @time tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],params)
        else
            tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],params)
        end

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k],ξ=ξ)
    end

    h5open(DetFilename(params,ξ=ξ), "w") do file
        # Active fraction   
        write(file,"xi",ξ)
        # Frequency grids
        write(file,"omega",real(ωlist))
        write(file,"eta",imag(ωlist))
        # Results
        write(file,"det",tabdetXi)
        # Parameters
        WriteParameters(file,params)
    end
    return tabdetXi
end



"""FindZeroCrossing(Ωguess,ηguess,FHT,params[,ξ,NITER,ACCURACY])

Newton-Raphson descent to find the zero crossing
"""
function FindZeroCrossing(Ωguess::Float64,ηguess::Float64,
                          FHT::FiniteHilbertTransform.AbstractFHT,
                          params::LinearParameters;
                          ξ::Float64=1.0,
                          NITER::Int64=32,
                          ACCURACY::Float64=1.0e-10)


    # Preparinng computations of the response matrices
    MMat, tabaMcoef, tabωminωmax = PrepareM(params)
    IMat = makeIMat(params.nradial)

    omgval = Ωguess + im*ηguess
    domega = 1.e-4

    completediterations = 0

    # this must be indicative of the multiprocessing bug: Threads helps here, even for 1
    #Threads.@threads for i = 1:NITER
    for i = 1:NITER

        # calculate the new off omega value
        omgvaloff = omgval + im*domega

        (params.VERBOSE > 1) && println("Step number $i: omega=$omgval, omegaoff=$omgvaloff")

        tabM!(omgval,MMat,tabaMcoef,tabωminωmax,FHT,params)

        centralvalue = detXi(IMat,MMat,ξ=ξ)

        tabM!(omgvaloff,MMat,tabaMcoef,tabωminωmax,FHT,params)

        offsetvalue = detXi(IMat,MMat,ξ=ξ)

        # ignore the imaginary part
        derivative = real(offsetvalue - centralvalue)/domega

        # take a step in omega given the derivative
        stepsize = real(centralvalue)/derivative
        omgval  = omgval - im*stepsize

        (params.VERBOSE > 1) && println("Newomg=$omgval, cval=$centralvalue, oval=$offsetvalue")

        # record iteration number
        completediterations += 1

        # check if we have gotten close enough already
        if abs(stepsize) < ACCURACY
            break
        end

    end

    (params.VERBOSE > 0) && println("LinearResponse.Xi.FindZeroCrossing: zero found in $completediterations steps.")

    return omgval
end