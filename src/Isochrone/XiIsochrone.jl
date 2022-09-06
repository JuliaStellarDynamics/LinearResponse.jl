"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

using LinearAlgebra


"""tabM!(omg,tabM,tabaMcoef,tabResVec,tabnpnq,structtabLeg,Ω0)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Ω0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabMIsochrone!(omg::Complex{Float64},
                        tabM::Array{Complex{Float64},2},
                        tabaMcoef::Array{Float64,4},
                        tabResVec::Matrix{Int64},
                        tabnpnq::Matrix{Int64},
                        FHT::FiniteHilbertTransform.FHTtype,
                        nradial::Int64,
                        Ω0::Float64=1.0)

    # get dimensions from the relevant tables
    nbnpnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]
    Ku      = FHT.Ku

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    #println("CallAResponse.Xi.tabM!: loop dimensions npnq=$nbnpnq, nResVec=$nbResVec, Ku=$Ku.")

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nResVec=1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec]

        # Rescale to get dimensionless frequency
        omgnodim = omg/Ω0

        # get the rescaled frequency
        ϖ = OrbitalElements.GetϖIsochrone(omgnodim,n1,n2)

        # get the Legendre integration values
        FiniteHilbertTransform.GettabD!(ϖ,FHT)

        # mame of the array where the Dk(w) are stored
        tabDLeg = structtabLeg.tabDLeg

        # loop over the basis indices to consider
        for inpnq=1:nbnpnq

            # get current value of (np,nq)
            np, nq = tabnpnq[1,inpnq], tabnpnq[2,inpnq]

            res = 0.0 + 0.0*im

            # loop over the Legendre functions to add all contributions
            for k=1:Ku
                res += tabaMcoef[nResVec,np,nq,k]*tabDLeg[k]
            end

            # fill the full M matrix:
            # as tabnpnq is the upper triangular matrix (with the diagonal),
            # we need to duplicate for symmetries

            # fill in the element (np,nq)
            tabM[np,nq] += res

            # fill in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.
            tabM[nq,np] += res
        end # basis index loop

        # do the actual writing
        #write(file, "tabMn1"*string(n1)*"n2"*string(n2),tabMtmp)

    end # resonance loop

    # contributions were added twice for diagonal elements: reset
    for np=1:nradial
        tabM[np,np] *= 0.5
    end


end




function RunMIsochrone(omglist::Array{Complex{Float64}},
                       gfuncdir::String,modedir::String,
                       Ku::Int64,Kv::Int64,Kw::Int64,
                       basis::AstroBasis.Basis_type,
                       FHT::FiniteHilbertTransform.FHTtype,
                       lharmonic::Int64,
                       n1max::Int64,
                       nradial::Int64,
                       Ω0::Float64,
                       modelname::String,dfname::String,
                       rb::Float64;
                       VERBOSE::Int64=0)

    nomglist = length(omglist)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(gfuncdir=gfuncdir,modedir=modedir)
    if checkdirs < 0
        return 0
    end

    # get basis parameters
    ndim = basis.dimension

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)


    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients ak
    MakeaMCoefficients(tabResVec,tabnpnq,FHT,gfuncdir,modedir,modelname,dfname,lharmonic,nradial,rb,VERBOSE=VERBOSE,OVERWRITE=false)

    # allocate structs for Dk(omega) computation
    structtabLeglist = [deepcopy(FHT) for k=1:Threads.nthreads()]

    # allocate memory for the response matrices M and identity matrices
    tabMlist = [zeros(Complex{Float64},nradial,nradial) for k=1:Threads.nthreads()]

    # make identity matrix and copies
    IMat = makeIMat(nradial)
    IMatlist = [deepcopy(IMat) for k=1:Threads.nthreads()]

    # allocate containers for determinant and min eigenvalue
    nomg = length(omglist)
    tabdetXi = zeros(Complex{Float64},nomg)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tabnpnq,Ku,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)
    println("CallAResponse.Xi.RunMIsochrone: tabaMcoef loaded.")

    println("CallAResponse.Xi.RunMIsochrone: computing $nomglist frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nomg

        k = Threads.threadid()

        if i==1
            @time tabMIsochrone!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tabnpnq,structtabLeglist[k],nradial,Ω0)
        else
            tabMIsochrone!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tabnpnq,structtabLeglist[k],nradial,Ω0)
        end

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])

    end

    WriteDeterminant(DetFilename(modedir,modelname,dfname,lharmonic,n1max,Ku,rb),omglist,tabdetXi)

    return tabdetXi
end



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTablesIsochrone(inputfile,
                                    omgval::Complex{Float64};
                                    VERBOSE::Int64=0)

    # include configuration parameters
    LoadConfiguration(inputfile)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(gfuncdir=gfuncdir,modedir=modedir)
    if checkdirs < 0
        return 0
    end

    # Construct the table of needed resonance vectors
    # Number of resonance vectors
    nbResVec = getnbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # get all Legendre weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(Ku)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients ak
    #MakeaMCoefficients(tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,rb=rb)
    MakeaMCoefficients(tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modedir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,OVERWRITE=false,rb=rb)

    # load aXi values
    tabaMcoef =  StageaMcoef(tabResVec,tabnpnq,Ku,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)
    println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # struct for Dk(omega) computation
    structtabLeglist = FiniteHilbertTransform.structtabLegcreate(Ku)

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1
    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    #tabM!(omgval,MMat,tabaMcoef,tabResVec,tabnpnq,structtabLeglist,dψ,d2ψ,nradial,Ω0)
    tabMIsochrone!(omgval,MMat,tabaMcoef,tabResVec,tabnpnq,structtabLeglist,nradial,Ω0)
    println("CallAResponse.Mode.ComputeModeTables: MMat constructed.")

    # eigenvalue, eigenfunction (eigenvector), eigenmode (for basis projection)
    EV,EF,EM = mevXi(MMat)

    return EV,EF,EM
end



"""
Newton-Raphson descent to find the zero crossing
"""
function FindZeroCrossingIsochrone(Ωguess::Float64,Etaguess::Float64,
                                   gfuncdir::String,modedir::String,
                                   Ku::Int64,Kv::Int64,Kw::Int64,
                                   basis::AstroBasis.Basis_type,
                                   lharmonic::Int64,
                                   n1max::Int64,
                                   nradial::Int64,
                                   Ω0::Float64,
                                   modelname::String,dfname::String,
                                   rb::Float64;
                                   NITER::Int64=32,
                                   eta::Bool=true,
                                   ACCURACY::Float64=1.0e-10,
                                   VERBOSE::Int64=0)

    # check directories
    if !(isdir(gfuncdir) && isdir(modedir))
        error("CallAResponse.Xi.jl: gfuncdir or modedir not found.")
    end

    # get basis parameters
    ndim = basis.dimension

    #####
    # Construct the table of needed resonance vectors
    #####
    nbResVec = getnbResVec(lharmonic,n1max,ndim) # Number of resonance vectors. ATTENTION, it is for the harmonics lharmonic
    tabResVec = maketabResVec(lharmonic,n1max,ndim) # Filling in the array of resonance vectors (n1,n2)

    # get all weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(Ku)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients ak
    MakeaMCoefficients(tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modedir,modelname,dfname,lharmonic,nradial,rb,VERBOSE=VERBOSE)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tabnpnq,Ku,nradial,
                                          modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)
    println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # Structs for Dk(omega) computation
    structtabLeglist = FiniteHilbertTransform.structtabLegcreate(Ku)
    # memory for the response matrices M and identity matrices
    tabMlist = zeros(Complex{Float64},nradial,nradial)
    IMat = makeIMat(nradial)

    omgval = Ωguess + im*Etaguess
    domega = 1.e-4

    completediterations = 0

    # this must be indicative of the multiprocessing bug: Threads helps here, even for 1
    #Threads.@threads for i = 1:NITER
    for i = 1:NITER

        # calculate the new off omega value
        omgvaloff = omgval + im*domega

        println("Step number $i: omega=$omgval, omegaoff=$omgvaloff")

        tabMIsochrone!(omgval,tabMlist,tabaMcoef,
              tabResVec,tabnpnq,
              structtabLeglist,
              nradial,Ω0,rmin,rmax)

        centralvalue = detXi(IMat,tabMlist)

        tabMIsochrone!(omgvaloff,tabMlist,tabaMcoef,
              tabResVec,tabnpnq,
              structtabLeglist,
              nradial,Ω0,rmin,rmax)

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
        println("CallAResponse.Xi.FindZeroCrossingIsochrone: zero found in $completediterations steps.")
    end

    return omgval

end




"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTablesIsochrone(omgval::Complex{Float64},
                                    gfuncdir::String,modedir::String,
                                    Ku::Int64,Kv::Int64,Kw::Int64,
                                    basis::AstroBasis.Basis_type,
                                    lharmonic::Int64,
                                    n1max::Int64,
                                    nradial::Int64,
                                    Ω0::Float64,
                                    modelname::String,dfname::String,
                                    rb::Float64,
                                    VERBOSE::Int64=0)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(gfuncdir=gfuncdir,modedir=modedir)
    if checkdirs < 0
        return 0
    end

    # Construct the table of needed resonance vectors
    # Number of resonance vectors
    nbResVec = getnbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # get all Legendre weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(Ku)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients ak
    MakeaMCoefficients(tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,rb,VERBOSE=VERBOSE)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tabnpnq,Ku,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)
    println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # struct for Dk(omega) computation
    structtabLeglist = FiniteHilbertTransform.structtabLegcreate(Ku)

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1
    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    tabMIsochrone!(omgval,MMat,tabaMcoef,tabResVec,tabnpnq,structtabLeglist,nradial,Omega0,rb)
    println("CallAResponse.Mode.ComputeModeTables: MMat constructed.")

    # eigenvalue, eigenfunction (eigenvector), eigenmode (for basis projection)
    EV,EF,EM = mevXi(MMat)

    return EV,EF,EM
end
