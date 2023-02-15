"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

using LinearAlgebra


"""makeaMCoefficients!(tabaMcoef,tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,basedir)

function to make the decomposition coefficients "a" of the response matrix M

these values do not depend on the frequency being evaluated: which makes them good to table

is this struggling from having to pass around a gigantic array? what if we did more splitting?
"""
function MakeaMCoefficients(FHT::FiniteHilbertTransform.FHTtype,
                            Parameters::ResponseParameters)

    # get relevant sizes
    nbResVec, tabResVec = Parameters.nbResVec, Parameters.tabResVec
    nradial  = Parameters.nradial
    Ku       = FHT.Ku

    (Parameters.VERBOSE > 2) && println("CallAResponse.Xi.MakeaMCoefficients: Check params: Ku=$Ku, nbResVec=$nbResVec.")

    # allocate the workspace
    tabaMcoef = zeros(Float64,Ku,nradial,nradial)

    # allocate quadrature result and warnflag tables
    restab  = zeros(Float64,Ku)
    warntab = zeros(Float64,Ku)


    # loop through all resonance vectors
    # make this parallel, but check that the loop is exited cleanly?
    # Threads.@threads for nres in 1:nbResVec
    for nres = 1:nbResVec

        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.MakeaMCoefficients: Starting on ($n1,$n2).")

        outputfilename = AxiFilename(n1,n2,Parameters)
        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,Parameters,"CallAResponse.Xi.MakeaMCoefficients: ($n1,$n2) resonance") || continue

        # Read G(u) from the GFunc file for this resonance
        gfuncfilename  = GFuncFilename(n1,n2,Parameters)
        tabGXi = h5read(gfuncfilename,"Gmat")

        for np = 1:nradial
            for nq = np:nradial
                # get the contribution
                FiniteHilbertTransform.GetaXi!(FHT,view(tabGXi,nq,np,:),restab,warntab)

                for k=1:Ku
                    # Warning if to many Inf or Nan values
                    (warntab[k] > 3) && println("CallAResponse.Xi.MakeaMCoefficients: NaN/Inf (warnflag=$(warntab[k])) values for (n1,n2)=($n1,$n2), (np,nq)=($np,$nq), and k=$k: $(restab[k]).")

                    # populate the symmetric matrix
                    tabaMcoef[k,nq,np] = restab[k] # Element (np,nq)
                    tabaMcoef[k,np,nq] = restab[k] # Element (nq,np). If np=nq overwrite (this is fine).
                end
            end
        end

        ωmin, ωmax = h5read(gfuncfilename,"omgmin"), h5read(gfuncfilename,"omgmax")
        # this is expensive enough to compute that we will want to save these
        # with the table fully constructed, loop back through to write after opening a file for the resonance
        h5open(outputfilename, "w") do outputfile
            # Mappings parameters
            write(outputfile, "omgmin",ωmin)
            write(outputfile, "omgmax",ωmax)
            # M_{n1,n2} decomposition coefficients
            write(outputfile,"aXi",tabaMcoef)
            # Parameters
            WriteParameters(outputfile,Parameters)
        end

    end # resonance vector loop

end



"""put all aXi values into memory"""
function StageaMcoef(Parameters::ResponseParameters)


    # get dimensions from the relevant tables
    nbResVec, tabResVec = Parameters.nbResVec, Parameters.tabResVec
    nradial  = Parameters.nradial
    Ku       = Parameters.Ku

    # allocate memory
    tabaMcoef = zeros(Float64,Ku,nradial,nradial,nbResVec)

    #Threads.@threads for nres=1:nbResVec # Loop over the resonances
    for nres = 1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres] # Current resonance (n1,n2)

        # retrieve the correct M table
        filename = AxiFilename(n1,n2,Parameters)
        tmptabaMcoef = h5read(filename,"aXi")

        for np = 1:nradial
            for nq = 1:nradial
                for k=1:Ku
                    # fill in the (np,nq) value
                    tabaMcoef[k,nq,np,nres] = tmptabaMcoef[k,nq,np]
                end
            end
        end

    end # end resonance vector loop

    return tabaMcoef
end



"""
    Stageωminωmax(Parameters::ResponseParameters)

    read the WMat HDF5 files to construct a matrix with the ωmin and ωmax values for every resonances
"""
function Stageωminωmax(Parameters::ResponseParameters)

    nbResVec, tabResVec = Parameters.nbResVec, Parameters.tabResVec
    tabωminωmax = zeros(Float64,2,nbResVec)
    
    for nres = 1:nbResVec
        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        # get already computed ωmin and ωmax values
        filename = AxiFilename(n1,n2,Parameters)
        ωmin = h5read(filename,"omgmin")
        ωmax = h5read(filename,"omgmax")

        tabωminωmax[1,nres], tabωminωmax[2,nres] = ωmin, ωmax
    end

    return tabωminωmax
end

"""tabM!(omg,tabM,tabaMcoef,FHT,Parameters)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Ω0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabM!(ω::Complex{Float64},
               tabM::Matrix{Complex{Float64}},
               tabaMcoef::Array{Float64,4},
               tabωminωmax::Matrix{Float64},
               FHT::FiniteHilbertTransform.FHTtype,
               Parameters::ResponseParameters)

    # get dimensions from the relevant tables
    nbResVec, tabResVec = Parameters.nbResVec, Parameters.tabResVec
    KuTruncation = Parameters.KuTruncation
    nradial  = Parameters.nradial
    VERBOSE  = Parameters.VERBOSE
    Ku       = FHT.Ku

    if KuTruncation < Ku
        Ku = KuTruncation
        (VERBOSE > 2) && println("CallAResponse.Xi.tabM!: truncating Ku series from $(FHT.Ku) to $Ku.")
    end

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    # Rescale to get dimensionless frequency
    ωnodim = ω/Parameters.OEparams.Ω₀

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nres = 1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        # get ωmin and ωmax values
        ωmin, ωmax = tabωminωmax[1,nres], tabωminωmax[2,nres]

        # get the rescaled frequency
        ϖ = OrbitalElements.Getϖ(ωnodim,ωmin,ωmax)

        # get the integration values
        FiniteHilbertTransform.GettabD!(ϖ,FHT)

        # mame of the array where the D_k(w) are stored
        tabD = FHT.tabDLeg

        # loop over the basis indices to consider
        for np = 1:nradial
            for nq = np:nradial

                res = 0.0 + 0.0*im

                # loop over the Legendre functions to add all contributions
                for k=1:Ku

                    # hard check for nans
                    val = tabaMcoef[k,nq,np,nres]*tabD[k]
                    if !isnan(val)
                        res += val
                    else
                        (k==1) && (VERBOSE>1) && println("CallAResponse.Xi.tabM!: NaN found for n=($n1,$n2), npnq=($np,$nq), k=$k")
                    end


                end

                # fill the full M matrix:
                # as tab_npnq is the upper triangular matrix (with the diagonal),
                # we need to duplicate for symmetries

                # fill in the element (np,nq)
                tabM[np,nq] += res

                # fill in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.
                tabM[nq,np] += res
            end
        end # basis index loop

    end # resonance loop

    # contributions were added twice for diagonal elements: reset
    for np=1:nradial
        tabM[np,np] *= 0.5
    end
end


"""
    detXi(IMat,tabM)

determinant of the susceptibility matrix I - M for known M.
"""
function detXi(IMat::Array{Complex{Float64},2},
               tabM::Array{Complex{Float64},2})

    # Computing the determinant of (I-M).
    # ATTENTION, we tell julia that the matrix is symmetric
    val = det(Symmetric(IMat-tabM))

    # only save the real portion
    return val # Output
end


"""
@TO DESCRIBE
"""
function RunM(ωlist::Array{Complex{Float64}},
              FHT::FiniteHilbertTransform.FHTtype,
              Parameters::ResponseParameters)

    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # check the basis values against the Parameters

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.RunM: Constructing M coefficients.")

    # make the decomposition coefficients a_k
    MakeaMCoefficients(FHT,Parameters)

    # allocate memory for FHT structs
    FHTlist = [deepcopy(FHT) for k=1:Threads.nthreads()]

    # allocate memory for the response matrices M and identity matrices
    tabMlist = [zeros(Complex{Float64},Parameters.nradial,Parameters.nradial) for k=1:Threads.nthreads()]

    # make identity matrix and copies
    IMat = makeIMat(Parameters.nradial)
    IMatlist = [deepcopy(IMat) for k=1:Threads.nthreads()]

    # how many omega values are we computing?
    nω = length(ωlist)
    # allocate containers for determinant and min eigenvalue
    tabdetXi = zeros(Complex{Float64},nω)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.RunM: Loading tabaMcoef...")

    # load aXi values
    tabaMcoef = StageaMcoef(Parameters)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.RunM: tabaMcoef loaded.")

    (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.RunM: computing $nω frequency values.")

    # load ωmin, ωmax values
    tabωminωmax = Stageωminωmax(Parameters)

    # loop through all frequencies
    Threads.@threads for i = 1:nω

        k = Threads.threadid()

        if (i==2) && (Parameters.VERBOSE>0) # skip the first in case there is compile time built in
            @time tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],Parameters)
        else
            tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabωminωmax,FHTlist[k],Parameters)
        end

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])
    end

    WriteDeterminant(DetFilename(Parameters),ωlist,tabdetXi)
    WriteParameters(DetFilename(Parameters),Parameters)
    return tabdetXi
end



"""FindZeroCrossing(Ωguess,ηguess,FHT,params[,NITER,eta,ACCURACY=1.0e-10])

Newton-Raphson descent to find the zero crossing
"""
function FindZeroCrossing(Ωguess::Float64,ηguess::Float64,
                          FHT::FiniteHilbertTransform.FHTtype,
                          Parameters::ResponseParameters;
                          NITER::Int64=32,
                          eta::Bool=true,
                          ACCURACY::Float64=1.0e-10)


    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(FHT,Parameters)

    # load aXi values
    tabaMcoef = StageaMcoef(Parameters)

    # load ωmin, ωmax values
    tabωminωmax = Stageωminωmax(Parameters)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # memory for the response matrices M and identity matrices
    tabMlist = zeros(Complex{Float64},Parameters.nradial,Parameters.nradial)
    IMat = makeIMat(Parameters.nradial)

    omgval = Ωguess + im*ηguess
    domega = 1.e-4

    completediterations = 0

    # this must be indicative of the multiprocessing bug: Threads helps here, even for 1
    #Threads.@threads for i = 1:NITER
    for i = 1:NITER

        # calculate the new off omega value
        omgvaloff = omgval + im*domega

        (Parameters.VERBOSE > 1) && println("Step number $i: omega=$omgval, omegaoff=$omgvaloff")

        tabM!(omgval,tabMlist,tabaMcoef,tabωminωmax,FHT,Parameters)

        centralvalue = detXi(IMat,tabMlist)

        tabM!(omgvaloff,tabMlist,tabaMcoef,tabωminωmax,FHT,Parameters)

        offsetvalue = detXi(IMat,tabMlist)

        # ignore the imaginary part
        derivative = real(offsetvalue - centralvalue)/domega

        # take a step in omega given the derivative
        stepsize = real(centralvalue)/derivative
        omgval  = omgval - im*stepsize

        (Parameters.VERBOSE > 1) && println("Newomg=$omgval, cval=$centralvalue, oval=$offsetvalue")

        # record iteration number
        completediterations += 1

        # check if we have gotten close enough already
        if abs(stepsize) < ACCURACY
            break
        end

    end

    (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.FindZeroCrossing: zero found in $completediterations steps.")

    return omgval
end