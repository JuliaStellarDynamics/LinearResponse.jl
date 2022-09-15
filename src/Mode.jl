



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTables(omgval::Complex{Float64},
                           dψ::Function,d2ψ::Function,
                           FHT::FiniteHilbertTransform.FHTtype,
                           gfuncdir::String,modedir::String,
                           Kv::Int64,
                           basis::AstroBasis.Basis_type,
                           lharmonic::Int64,
                           n1max::Int64,
                           Ω₀::Float64,
                           rmin::Float64,rmax::Float64,
                           modelname::String,dfname::String;
                           VERBOSE::Int64=0)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(gfuncdir=gfuncdir,modedir=modedir)
    if checkdirs < 0
        return 0
    end

    # get needed parameters from structures
    Ku      = FHT.Ku
    nradial = basis.nmax
    rb      = basis.rb
    ndim    = basis.dimension

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tabnpnq,FHT,gfuncdir,modedir,modelname,dfname,lharmonic,nradial,rb,Kv,VERBOSE=VERBOSE)

    # load aXi values
    tabaMcoef = StageaMcoef(tabResVec,tabnpnq,FHT.Ku,Kv,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)

    if VERBOSE>0
        println("CallAResponse.Xi.ComputeModeTables: tabaMcoef loaded.")
    end

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1
    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    tabM!(omgval,MMat,tabaMcoef,tabResVec,tabnpnq,FHT,dψ,d2ψ,nradial,Ω₀,rmin,rmax,VERBOSE=VERBOSE)

    if VERBOSE>0
        println("CallAResponse.Mode.ComputeModeTables: MMat constructed.")
    end

    # eigenvalue, eigenfunction (eigenvector), eigenmode (for basis projection)
    EV,EF,EM = mevXi(MMat)

    return EV,EF,EM

end


"""
    mevXi(tabM)

minimal eigenvalue of M
"""
function mevXi(tabM::Array{Complex{Float64},2})

    # these should be equal, and =nradial!
    nEig1,nEig2 = size(tabM)

    # Computing the eigenvalue that is closest to 1
    tabeigvals = eigvals(tabM)
    tabeigvecs = eigvecs(tabM)

    indEig = 1
    for indrad = 1:nEig1

        # did we find an eigenvalue that is even closer to 1.0?
        if (abs(1.0-tabeigvals[indrad]) < abs(1.0-tabeigvals[indEig]))
            # if yes, update the index of the eigenvalue
            indEig = indrad
        end
    end

    # construct the mode table
    tabEigenMode = zeros(Float64,nEig1)

    # now fill in the eigenmode
    # loop over the number of basis elements
    for np=1:nEig1

        # extract the eigenvector
        # tabeigvecs is the matrix whose COLUMNS are eigenvectors
        # careful, we are only getting the real part
        tabEigenMode[np] = real(tabeigvecs[np,indEig])
    end

    # output
    return tabeigvals[indEig],tabeigvecs[:,indEig],tabEigenMode

end




"""
    GetModeShape(basis,l,Rmin,Rmax,nR,EigenMode)
Function that computes the radial shape of a given mode
-see Eq. (72) of Hamilton+ (2018)



"""
function GetModeShape(basis::AstroBasis.Basis_type,
                      lharmonic::Int64,
                      nradial::Int64,
                      n1max::Int64,
                      Rmin::Float64,Rmax::Float64,
                      nRMode::Int64,
                      EigenMode::Array{Float64},
                      modedir::String,
                      modelname::String,
                      dfname::String,
                      Ku::Int64;
                      VERBOSE::Int64=0)

    # prep the basis
    AstroBasis.fill_prefactors!(basis)

    # get the step distance of the array in RMode
    deltaRMode = (Rmax - Rmin)/(nRMode - 1)

    # table of R for which the mode is computed
    radiusvals = LinRange(Rmin,Rmax,nRMode)

    # table containing the radial shape of the mode
    ModeRadius         = zeros(Float64,nRMode)
    ModePotentialShape = zeros(Float64,nRMode)
    ModeDensityShape   = zeros(Float64,nRMode)

    println("CallAResponse.Mode.GetModeShape: Starting radius loop...")

    # at each radius, compute the shape of the mode response
    for irad=1:nRMode

        # current value of R
        R = radiusvals[irad]

        # initialise the values
        pval = 0.0
        dval = 0.0

        # for each basis element, add the coefficient contribution
        for np=1:nradial

            # add the contribution from the basis elements.
            pval += EigenMode[np]*AstroBasis.getUln(basis,lharmonic,np-1,R)
            dval += EigenMode[np]*AstroBasis.getDln(basis,lharmonic,np-1,R)

        end

        # log contribution in tables
        ModeRadius[irad] = R
        ModePotentialShape[irad] = pval
        ModeDensityShape[irad] = dval

    end

    h5open(ModeFilename(modedir,modelname,dfname,lharmonic,n1max,Ku), "w") do file
        write(file,"ModeRadius",ModeRadius)         # write tabRMode to file
        write(file,"ModePotentialShape",ModePotentialShape) # write tabShapeMode to file
        write(file,"ModeDensityShape",ModeDensityShape) # write tabShapeMode to file
    end

    # return just in case we want to do something else
    return ModeRadius,ModePotentialShape,ModeDensityShape

end

"""
    GetModeShapeComplex(basis,l,Rmin,Rmax,nR,EigenModeC)
Function that computes the radial and azimuthal shape of a given mode
-see Eq. (72) of Hamilton+ (2018)
"""
function GetModeShapeComplex(basis::AstroBasis.Basis_type,
                            lharmonic::Int64,
                            n1max::Int64,
                            Rmin::Float64,Rmax::Float64,
                            nRMode::Int64,
                            EigenMode::Array{Complex{Float64}},
                            modedir::String,
                            modelname::String,
                            dfname::String,
                            Ku::Int64;
                            VERBOSE::Int64=0)

    # prep the basis
    AstroBasis.fill_prefactors!(basis)

    # Basis parameters
    nradial = basis.nmax

    # get the step distance of the array in RMode
    deltaRMode = (Rmax - Rmin)/(nRMode - 1)

    # table of R for which the mode is computed
    radiusvals = LinRange(Rmin,Rmax,nRMode)

    # table containing the radial shape of the mode
    ModeRadius         = zeros(Float64,nRMode)
    ModePotentialShape = zeros(Complex{Float64},nRMode)
    ModeDensityShape   = zeros(Complex{Float64},nRMode)

    println("CallAResponse.Mode.GetModeShape: Starting radius loop...")

    # at each radius, compute the shape of the mode response
    for irad=1:nRMode

        # current value of R
        R = radiusvals[irad]

        # initialise the values
        pval = 0.0
        dval = 0.0

        # for each basis element, add the coefficient contribution
        for np=1:nradial

            # add the contribution from the basis elements.
            pval += EigenMode[np]*AstroBasis.getUln(basis,lharmonic,np-1,R)
            dval += EigenMode[np]*AstroBasis.getDln(basis,lharmonic,np-1,R)

        end

        # log contribution in tables
        ModeRadius[irad] = R
        ModePotentialShape[irad] = pval
        ModeDensityShape[irad] = dval

    end

    h5open(ModeFilename(modedir,modelname,dfname,lharmonic,n1max,Ku), "w") do file
        write(file,"ModeRadius",ModeRadius)         # write tabRMode to file
        write(file,"ModePotentialShape",ModePotentialShape) # write tabShapeMode to file
        write(file,"ModeDensityShape",ModeDensityShape) # write tabShapeMode to file
    end

    # return just in case we want to do something else
    return ModeRadius,ModePotentialShape,ModeDensityShape

end
