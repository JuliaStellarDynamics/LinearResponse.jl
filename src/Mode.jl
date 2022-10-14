



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTables(omgval::Complex{Float64},
                           FHT::FiniteHilbertTransform.FHTtype,
                           basis::AstroBasis.Basis_type,
                           Parameters::ResponseParameters)

    # Check directory names
    CheckConfigurationDirectories([Parameters.gfuncdir,Parameters.modedir]) || (return 0)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(FHT,Parameters)

    # load aXi values
    tabaMcoef = StageaMcoef(Parameters)

    if Parameters.VERBOSE>0
        println("CallAResponse.Xi.ComputeModeTables: tabaMcoef loaded.")
    end

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},Parameters.nradial,Parameters.nradial)

    tabM!(omgval,MMat,tabaMcoef,FHT,Parameters)

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
function mevXi(tabM::Matrix{Complex{Float64}})

    # these should be equal, and =nradial!
    nEig1, _ = size(tabM)

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
    tabEigenMode = zeros(Complex{Float64},nEig1)

    # now fill in the eigenmode
    # loop over the number of basis elements
    for np=1:nEig1

        # extract the eigenvector
        # tabeigvecs is the matrix whose COLUMNS are eigenvectors
        tabEigenMode[np] = tabeigvecs[np,indEig]
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
                      EigenMode::Vector,
                      Parameters::ResponseParameters)

    # table of R for which the mode is computed
    radiusvals = LinRange(Rmin,Rmax,nRMode)

    # table containing the radial shape of the mode
    eigentype = typeof(EigenMode[1])
    ModeRadius         = zeros(Float64,nRMode)
    ModePotentialShape = zeros(eigentype,nRMode)
    ModeDensityShape   = zeros(eigentype,nRMode)

    (Parameters.VEBOSE > 1) && println("CallAResponse.Mode.GetModeShape: Starting radius loop...")

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