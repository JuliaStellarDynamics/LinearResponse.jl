



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTables(omgval::Complex{Float64},
                           dψ::Function,d2ψ::Function,
                           FHT::FiniteHilbertTransform.FHTtype,
                           basis::AstroBasis.Basis_type,
                           Parameters::ResponseParameters)

    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # get needed parameters from structures
    Ku      = FHT.Ku
    nradial = basis.nmax
    rb      = basis.rb
    ndim    = basis.dimension

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(Parameters.lharmonic,Parameters.n1max,ndim)


    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,FHT,Parameters)

    # load aXi values
    tabaMcoef = StageaMcoef(tabResVec,Parameters)

    if Parameters.VERBOSE>0
        println("CallAResponse.Xi.ComputeModeTables: tabaMcoef loaded.")
    end

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1
    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    tabM!(omgval,MMat,tabaMcoef,tabResVec,FHT,dψ,d2ψ,nradial,Parameters.Ω₀,Parameters.rmin,Parameters.rmax,VERBOSE=Parameters.VERBOSE)

    if Parameters.VERBOSE>0
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
                      Rmin::Float64,Rmax::Float64,
                      nRMode::Int64,
                      EigenMode::Array{Float64},
                      Parameters::ResponseParameters)

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
        for np=1:Parameters.nradial

            # add the contribution from the basis elements.
            pval += EigenMode[np]*AstroBasis.getUln(basis,Parameters.lharmonic,np-1,R)
            dval += EigenMode[np]*AstroBasis.getDln(basis,Parameters.lharmonic,np-1,R)

        end

        # log contribution in tables
        ModeRadius[irad] = R
        ModePotentialShape[irad] = pval
        ModeDensityShape[irad] = dval

    end

    h5open(ModeFilename(Parameters), "w") do file
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
                            Rmin::Float64,Rmax::Float64,
                            nRMode::Int64,
                            EigenMode::Array{Complex{Float64}},
                            Parameters::ResponseParameters)

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
            pval += EigenMode[np]*AstroBasis.getUln(basis,Parameters.lharmonic,np-1,R)
            dval += EigenMode[np]*AstroBasis.getDln(basis,Parameters.lharmonic,np-1,R)

        end

        # log contribution in tables
        ModeRadius[irad] = R
        ModePotentialShape[irad] = pval
        ModeDensityShape[irad] = dval

    end

    h5open(ModeFilename(Parameters), "w") do file
        write(file,"ModeRadius",ModeRadius)         # write tabRMode to file
        write(file,"ModePotentialShape",ModePotentialShape) # write tabShapeMode to file
        write(file,"ModeDensityShape",ModeDensityShape) # write tabShapeMode to file
    end

    # return just in case we want to do something else
    return ModeRadius,ModePotentialShape,ModeDensityShape

end
