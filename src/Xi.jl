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
                            tabwGLquad::Vector{Float64},
                            tabPGLquad::Matrix{Float64},
                            tabINVcGLquad::Vector{Float64},
                            gfuncdir::String,
                            modelname::String,
                            dfname::String,
                            lharmonic::Int64,
                            nradial::Int64;
                            VERBOSE::Int64=0,
                            OVERWRITE::Bool=false,
                            modedir::String="")

    # get relevant sizes
    K_u      = size(tabwGLquad)[1]
    nb_npnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    if VERBOSE > 2
        println("CallAResponse.Xi.MakeMCoefficients: Check params: K_u=$K_u, nb_npnq=$nb_npnq, nbResVec=$nbResVec.")
    end

    # allocate the workspace
    tabaMcoef = zeros(Float64,nradial,nradial,K_u)

    # loop through all resonance vectors
    # make this parallel, but check that the loop is exited cleanly?
    Threads.@threads for nresvec in 1:nbResVec
    #for nresvec in 1:nbResVec

        n1,n2 = tabResVec[1,nresvec],tabResVec[2,nresvec]

        # don't do this loop if the file calculation already exists (unless asked)
        outputfilename = AxiFilename(modedir,modelname,dfname,lharmonic,n1,n2,K_u,rb)
        if isfile(outputfilename)

            # log if requested
            if VERBOSE>0
                println("CallAResponse.Xi.MakeMCoefficients: file already exists for step $nresvec of $nbResVec, ($n1,$n2).")
            end

            # decide if we want to overwrite anyway
            if OVERWRITE
                println("...recomputing anyway.")
            else
                continue
            end

        end

        if VERBOSE > 0
            println("CallAResponse.Xi.MakeMCoefficients: on step $nresvec of $nbResVec: ($n1,$n2).")
        end

        # open the resonance file
        filename = GFuncFilename(gfuncdir,modelname,dfname,lharmonic,n1,n2,K_u,rb)
        inputfile = h5open(filename,"r")

        if VERBOSE > 0
            println("CallAResponse.Xi.MakeMCoefficients: opened file $filename.")
        end

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)

            # open the correct resonance vector
            # read in the correct G(u) function
            tabGXi = read(inputfile,"GXinp"*string(np)*"nq"*string(nq))

            # Loop over the Legendre functions
            for k=1:K_u

                res = 0.0 # Initialisation of the result

                for i=1:K_u # Loop over the G-L nodes

                    G = tabGXi[i] # Current value of G[u_i]

                    # check for NaN contribution: skip this contribution in that case
                    if isnan(G)
                        if VERBOSE > 1
                            println("CallAResponse.Xi.MakeMCoefficients: NaN value for (n1,n2,np,nq)=($n1,$n2,$np,$nq) and K_u=$i (of $K_u).")
                        end
                        continue
                    end

                    # check for INF contribution: skip the contribution in that case
                    if isinf(G)
                        if VERBOSE > 1
                            println("CallAResponse.Xi.MakeMCoefficients: Inf value for (n1,n2,np,nq)=($n1,$n2,$np,$nq) and K_u=$i (of $K_u).")
                        end
                        continue
                    end

                    w = tabwGLquad[i] # Current weight
                    P = tabPGLquad[k,i] # Current value of P_k. ATTENTION, to the order of the arguments.
                    res += w*G*P # Update of the sum
                end

                res *= tabINVcGLquad[k] # Multiplying by the Legendre prefactor.

                # populate the symmetric matrix
                tabaMcoef[np,nq,k] = res # Element (np,nq)
                tabaMcoef[nq,np,k] = res # Element (nq,np). If np=nq overwrite (this is fine).

            end

        end # basis function loop

        # close the Gfunc file
        close(inputfile)

        # this is expensive enough to compute that we will want to save these
        # with the table fully constructed, loop back through to write after opening a file for the resonance
        h5open(outputfilename, "w") do outputfile
            for i_npnq=1:nb_npnq
                np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)
                write(outputfile, "aXinp"*string(np)*"nq"*string(nq),tabaMcoef[np,nq,:])
            end
        end

    end # resonance vector loop

end



"""put all aXi values into memory"""
function StageaMcoef(tabResVec::Matrix{Int64},
                     tab_npnq::Matrix{Int64},
                     K_u::Int64,
                     nradial::Int64;
                     modedir::String="",
                     modelname::String="",
                     dfname::String="",
                     lharmonic::Int64="")


    # get dimensions from the relevant tables
    nb_npnq  = size(tab_npnq)[2]
    nbResVec = size(tabResVec)[2]

    # allocate memory
    tabaMcoef = zeros(Float64,nbResVec,nradial,nradial,K_u)


    #Threads.@threads for nResVec=1:nbResVec # Loop over the resonances
    for nResVec=1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec] # Current resonance (n1,n2)

        # retrieve the correct M table
        filename = AxiFilename(modedir,modelname,dfname,lharmonic,n1,n2,K_u,rb)
        inputfile = h5open(filename,"r")

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tab_npnq[1,i_npnq], tab_npnq[2,i_npnq] # Current value of (np,nq)

            # get the table from the inputfile
            tmptabaMcoef = read(inputfile,"aXinp"*string(np)*"nq"*string(nq))

            for k=1:K_u

                # fill in the (np,nq) value
                tabaMcoef[nResVec,np,nq,k] = tmptabaMcoef[k]

                # fill the other symmetric side of the array:
                # does not matter because we symmetrise later
                tabaMcoef[nResVec,nq,np,k] = tmptabaMcoef[k]

            end

        end # end basis loop

    end # end resonance vector loop

    return tabaMcoef


end

"""tabM!(omg,tabM,tabaMcoef,tabResVec,tabnpnq,struct_tabLeg,dψ,d2ψ,LINEAR,Ω0)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Ω0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabM!(omg::Complex{Float64},
               tabM::Array{Complex{Float64},2},
               tabaMcoef::Array{Float64,4},
               tabResVec::Matrix{Int64},
               tab_npnq::Matrix{Int64},
               struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
               dψ::Function,
               d2ψ::Function,
               nradial::Int64,
               LINEAR::String="unstable",
               Ω0::Float64=1.0)

    # get dimensions from the relevant tables
    nb_npnq  = size(tab_npnq)[2]
    nbResVec = size(tabResVec)[2]
    K_u      = size(tabaMcoef)[4]

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    #println("CallAResponse.Xi.tabM!: loop dimensions npnq=$nb_npnq, nResVec=$nbResVec, K_u=$K_u.")

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nResVec=1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec]

        # Rescale to get dimensionless frequency
        omg_nodim = omg/Ω0

        # get the rescaled frequency
        varpi = OrbitalElements.GetVarpi(omg_nodim,n1,n2,dψ,d2ψ,Ω₀=Ω0)

        # get the Legendre integration values
        PerturbPlasma.get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR)

        # mame of the array where the D_k(w) are stored
        tabDLeg = struct_tabLeg.tabDLeg

        # loop over the basis indices to consider
        for i_npnq=1:nb_npnq

            # get current value of (np,nq)
            np, nq = tab_npnq[1,i_npnq], tab_npnq[2,i_npnq]

            res = 0.0 + 0.0*im

            # loop over the Legendre functions to add all contributions
            for k=1:K_u

                # hard check for nans
                val = tabaMcoef[nResVec,np,nq,k]*tabDLeg[k]
                if !isnan(val)
                    res += val
                else
                    if k==1
                        println("CallAResponse.Xi.tabM!: NaN found for n=($n1,$n2), npnq=($np,$nq), k=$k")
                    end
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



function RunM(inputfile::String,
              omglist::Array{Complex{Float64}};
              VERBOSE::Int64=0)

    # need some sort of 'if' for whether this already exists
    include(inputfile)

    nomglist = length(omglist)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(gfuncdir=gfuncdir,modedir=modedir)
    if checkdirs < 0
        return 0
    end

    # calculate the number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # get all weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = PerturbPlasma.tabGLquad(K_u)

    # make the (np,nq) vectors that we need to evaluate
    tab_npnq = makeTabnpnq(nradial)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,modedir=modedir)

    # allocate structs for D_k(omega) computation
    struct_tabLeglist = [PerturbPlasma.struct_tabLeg_create(K_u) for k=1:Threads.nthreads()]

    # allocate memory for the response matrices M and identity matrices
    tabMlist = [zeros(Complex{Float64},nradial,nradial) for k=1:Threads.nthreads()]

    # make identity matrix and copies
    IMat = makeIMat(nradial)
    IMatlist = [deepcopy(IMat) for k=1:Threads.nthreads()]

    # allocate containers for determinant and min eigenvalue
    nomg = length(omglist)
    tabdetXi = zeros(Complex{Float64},nomg)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic)
    println("CallAResponse.Xi.RunM: tabaMcoef loaded.")

    println("CallAResponse.Xi.RunM: Starting frequency analysis, using $LINEAR integration.")
    println("CallAResponse.Xi.RunM: computing $nomglist frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nomg

        k = Threads.threadid()

        if i==2 # skip the first in case there is compile time built in
            @time tabM!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist[k],dψ,d2ψ,nradial,LINEAR,Ω0)
        else
            tabM!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist[k],dψ,d2ψ,nradial,LINEAR,Ω0)
        end

        # do we need some sort of diagnostic check here?
        #println("Mval=")

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])

    end

    WriteDeterminant(DetFilename(modedir,modelname,dfname,lharmonic,n1max,K_u,rb),omglist,tabdetXi)

    return tabdetXi
end



"""
Newton-Raphson descent to find the zero crossing
"""
function FindZeroCrossing(inputfile::String,
                          Ωguess::Float64,
                          Etaguess::Float64;
                          NITER::Int64=32,
                          eta::Bool=true,
                          Ω0::Float64=1.0,
                          ACCURACY::Float64=1.0e-10,
                          VERBOSE::Int64=0)

    # need some sort of 'if' for whether this already exists
    LoadConfiguration(inputfile)

    #####
    # Check directories names
    #####
    if !(isdir(gfuncdir) && isdir(modedir))
        error("CallAResponse.Xi.jl: gfuncdir or modedir not found.")
    end

    #####
    # Construct the table of needed resonance vectors
    #####
    nbResVec = get_nbResVec(lharmonic,n1max,ndim) # Number of resonance vectors. ATTENTION, it is for the harmonics lharmonic
    tabResVec = maketabResVec(lharmonic,n1max,ndim) # Filling in the array of resonance vectors (n1,n2)

    # get all weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = PerturbPlasma.tabGLquad(K_u)

    # make the (np,nq) vectors that we need to evaluate
    tab_npnq = makeTabnpnq(nradial)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,modedir=modedir)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial,
                                          modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic)
    println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # Structs for D_k(omega) computation
    struct_tabLeglist = PerturbPlasma.struct_tabLeg_create(K_u)
    # memory for the response matrices M and identity matrices
    tabMlist = zeros(Complex{Float64},nradial,nradial)
    IMat = makeIMat(nradial)

    omgval = Ωguess + im*Etaguess
    domega = 1.e-4

    completediterations = 0
    for i = 1:NITER

        # calculate the new off omega value
        omgvaloff = omgval + im*domega

        println("Step number $i: omega=$omgval, omegaoff=$omgvaloff")

        tabM!(omgval,tabMlist,tabaMcoef,
              tabResVec,tab_npnq,
              struct_tabLeglist,
              dψ,d2ψ,nradial,LINEAR,Ω0)

        centralvalue = detXi(IMat,tabMlist)

        tabM!(omgvaloff,tabMlist,tabaMcoef,
              tabResVec,tab_npnq,
              struct_tabLeglist,
              dψ,d2ψ,nradial,LINEAR,Ω0)

        offsetvalue = detXi(IMat,tabMlist)

        # ignore the imaginary part
        derivative = real(offsetvalue - centralvalue)/domega

        # take a step in omega given the derivative
        stepsize = real(centralvalue)/derivative
        omgval  = omgval - im*stepsize

        println("Newomg=$omgval, cval=$centralvalue, oval=$offsetvalue")

        # record iteration number
        completediterations += 1

        # check if we have gotten close enough already
        if abs(stepsize) < ACCURACY
            break
        end

    end

    if VERBOSE > 0
        println("CallAResponse.Xi.FindZeroCrossing: zero found in $completediterations steps.")
    end

    return omgval

end
