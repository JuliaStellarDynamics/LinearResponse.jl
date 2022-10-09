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
function MakeaMCoefficients(tabResVec::Matrix{Int64},
                            tabnpnq::Matrix{Int64},
                            FHT::FiniteHilbertTransform.FHTtype,
                            Parameters::ResponseParameters)

    # get relevant sizes
    Ku       = FHT.Ku
    nbnpnq   = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    (Parameters.VERBOSE > 2) && println("CallAResponse.Xi.MakeaMCoefficients: Check params: K_u=$K_u, nb_npnq=$nb_npnq, nbResVec=$nbResVec.")

    # allocate the workspace
    tabaMcoef = zeros(Float64,Parameters.nradial,Parameters.nradial,Ku)

    # loop through all resonance vectors
    # make this parallel, but check that the loop is exited cleanly?
    # Threads.@threads for nres in 1:nbResVec
    for nres in 1:nbResVec

        n1,n2 = tabResVec[1,nres],tabResVec[2,nres]

        # don't do this loop if the file calculation already exists (unless asked)
        outputfilename = AxiFilename(n1,n2,Parameters)
        if isfile(outputfilename)

            # log if requested
            (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.MakeaMCoefficients: file already exists for step $nres of $nbResVec, ($n1,$n2).")

            # decide if we want to overwrite anyway
            Parameters.OVERWRITE ? println("...recomputing anyway.") : continue
        end

        (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.MakeaMCoefficients: on step $nres of $nbResVec: ($n1,$n2).")

        # open the resonance file
        filename = GFuncFilename(n1,n2,Parameters)
        inputfile = h5open(filename,"r")
        tabGXi = read(inputfile,"GXinp")
        close(inputfile)

        (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.MakeaMCoefficients: opened file $filename.")

        # Loop over the basis indices to consider
        for i_npnq=1:nbnpnq
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)


            # get the contribution
            res,warnflag = FiniteHilbertTransform.GetaXi(FHT,tabGXi[np,nq,:])

            for k=1:Ku

                if warnflag[k] > 2
                    println("CallAResponse.Xi.MakeaMCoefficients: NaN/Inf (warnflag=$(warnflag[k])) values for (n1,n2)=($n1,$n2), (np,nq)=($np,$nq), and k=$k: $(res[k]).")
                end

                # populate the symmetric matrix
                tabaMcoef[np,nq,k] = res[k] # Element (np,nq)
                tabaMcoef[nq,np,k] = res[k] # Element (nq,np). If np=nq overwrite (this is fine).

            end

        end # basis function loop



        # this is expensive enough to compute that we will want to save these
        # with the table fully constructed, loop back through to write after opening a file for the resonance
        h5open(outputfilename, "w") do outputfile
            for i_npnq=1:nbnpnq
                np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)
                write(outputfile, "aXinp"*string(np)*"nq"*string(nq),tabaMcoef[np,nq,:])
            end
        end

    end # resonance vector loop

end



"""put all aXi values into memory"""
function StageaMcoef(tabResVec::Matrix{Int64},
                     tabnpnq::Matrix{Int64},
                     Parameters::ResponseParameters)


    # get dimensions from the relevant tables
    nbnpnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    # allocate memory
    tabaMcoef = zeros(Float64,nbResVec,Parameters.nradial,Parameters.nradial,Parameters.Ku)


    #Threads.@threads for nres=1:nbResVec # Loop over the resonances
    for nres=1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres] # Current resonance (n1,n2)

        # retrieve the correct M table
        filename = AxiFilename(n1,n2,Parameters)
        inputfile = h5open(filename,"r")

        # Loop over the basis indices to consider
        for i_npnq=1:nbnpnq
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)

            # get the table from the inputfile
            tmptabaMcoef = read(inputfile,"aXinp"*string(np)*"nq"*string(nq))

            for k=1:Parameters.Ku

                # fill in the (np,nq) value
                tabaMcoef[nres,np,nq,k] = tmptabaMcoef[k]

                # fill the other symmetric side of the array:
                # does not matter because we symmetrise later
                tabaMcoef[nres,nq,np,k] = tmptabaMcoef[k]

            end

        end # end basis loop

    end # end resonance vector loop

    return tabaMcoef


end


"""tabM!(omg,tabM,tabaMcoef,tabResVec,tabnpnq,struct_tabLeg,dψ,d2ψ,Ω0)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Ω0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabM!(ω::Complex{Float64},
               tabM::Array{Complex{Float64},2},
               tabaMcoef::Array{Float64,4},
               tabResVec::Matrix{Int64},
               tabnpnq::Matrix{Int64},
               FHT::FiniteHilbertTransform.FHTtype,
               dψ::Function,
               d2ψ::Function,
               nradial::Int64,
               Ω₀::Float64,
               rmin::Float64,rmax::Float64;
               VERBOSE::Int64=0,
               KuTruncation::Int64=10000)

    # get dimensions from the relevant tables
    nbnpnq      = size(tabnpnq)[2]
    nbResVec    = size(tabResVec)[2]
    Ku          = FHT.Ku

    if KuTruncation < Ku
        Ku = KuTruncation
        (VERBOSE > 2) && println("CallAResponse.Xi.tabM!: truncating Ku series from $(FHT.Ku) to $Ku.")
    end

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    # Rescale to get dimensionless frequency
    ωnodim = ω/Ω₀

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nres=1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        # get the rescaled frequency
        ϖ = OrbitalElements.Getϖ(ωnodim,n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

        # get the integration values
        FiniteHilbertTransform.GettabD!(ϖ,FHT)

        # mame of the array where the D_k(w) are stored
        tabD = FHT.tabDLeg


        # loop over the basis indices to consider
        for i_npnq=1:nbnpnq

            # get current value of (np,nq)
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq]

            res = 0.0 + 0.0*im

            # loop over the Legendre functions to add all contributions
            for k=1:Ku

                # hard check for nans
                val = tabaMcoef[nres,np,nq,k]*tabD[k]
                
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



function RunM(ωlist::Array{Complex{Float64}},
              dψ::Function,d2ψ::Function,
              FHT::FiniteHilbertTransform.FHTtype,
              basis::AstroBasis.Basis_type,
              Parameters::ResponseParameters
              )

    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # get basis parameters
    ndim, nradial, rb = basis.dimension, basis.nmax, basis.rb

    #if Parameters.nradial != nradial, break

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(Parameters.lharmonic,Parameters.n1max,ndim)

    # get all weights
    tabu, Ku = FHT.tabu, FHT.Ku

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(Parameters.nradial)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.RunM: Constructing M coefficients.")

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tabnpnq,FHT,Parameters)

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
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tabnpnq,Parameters)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.RunM: tabaMcoef loaded.")

    (Parameters.VERBOSE > 0) && println("CallAResponse.Xi.RunM: computing $nω frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nω

        k = Threads.threadid()

        if (i==2) && (Parameters.VERBOSE>0) # skip the first in case there is compile time built in
            @time tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabResVec,tabnpnq,FHTlist[k],dψ,d2ψ,
                        Parameters.nradial,Parameters.Ω₀,
                        Parameters.rmin,Parameters.rmax,
                        VERBOSE=Parameters.VERBOSE,KuTruncation=Parameters.KuTruncation)
        else
            tabM!(ωlist[i],tabMlist[k],tabaMcoef,tabResVec,tabnpnq,FHTlist[k],dψ,d2ψ,
                  Parameters.nradial,Parameters.Ω₀,
                  Parameters.rmin,Parameters.rmax,
                  VERBOSE=Parameters.VERBOSE,KuTruncation=Parameters.KuTruncation)
        end

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])

    end

    WriteDeterminant(DetFilename(Parameters),ωlist,tabdetXi)

    return tabdetXi
end



"""FindZeroCrossing(Ωguess,ηguess,dψ,d2ψ,gfuncdir,modedir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,dfname,rb,rmin::Float64,rmax[,NITER,eta,ACCURACY=1.0e-10,VERBOSE=0,OVERWRITE=false])

Newton-Raphson descent to find the zero crossing
"""
function FindZeroCrossing(Ωguess::Float64,ηguess::Float64,
                          dψ::Function,d2ψ::Function,
                          FHT::FiniteHilbertTransform.FHTtype,
                          basis::AstroBasis.Basis_type,
                          Parameters::ResponseParameters;
                          NITER::Int64=32,
                          eta::Bool=true,
                          ACCURACY::Float64=1.0e-10)


    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # get basis parameters
    ndim, nradial = basis.dimension, basis.nmax

    #####
    # Construct the table of needed resonance vectors
    #####
    nbResVec, tabResVec = MakeTabResVec(Parameters.lharmonic,Parameters.n1max,ndim) # Number of resonance vectors. ATTENTION, it is for the harmonics lharmonic

    # get all weights
    tabu, Ku = FHT.tabu, FHT.Ku

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tabnpnq,FHT,Parameters)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,Parameters)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # memory for the response matrices M and identity matrices
    tabMlist = zeros(Complex{Float64},nradial,nradial)
    IMat = makeIMat(nradial)

    omgval = Ωguess + im*ηguess
    domega = 1.e-4

    completediterations = 0

    # this must be indicative of the multiprocessing bug: Threads helps here, even for 1
    #Threads.@threads for i = 1:NITER
    for i = 1:NITER

        # calculate the new off omega value
        omgvaloff = omgval + im*domega

        (Parameters.VERBOSE > 1) && println("Step number $i: omega=$omgval, omegaoff=$omgvaloff")

        tabM!(omgval,tabMlist,tabaMcoef,
              tabResVec,tabnpnq,
              FHT,
              dψ,d2ψ,Parameters.nradial,Parameters.Ω₀,Parameters.rmin,Parameters.rmax)

        centralvalue = detXi(IMat,tabMlist)

        tabM!(omgvaloff,tabMlist,tabaMcoef,
              tabResVec,tab_npnq,
              FHT,
              dψ,d2ψ,Parameters.nradial,Parameters.Ω₀,Parameters.rmin,Parameters.rmax)

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
