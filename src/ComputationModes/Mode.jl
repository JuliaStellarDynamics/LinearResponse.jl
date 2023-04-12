



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTables(omgval::ComplexF64,
                           FHT::FiniteHilbertTransform.AbstractFHT,
                           params::LinearParameters)

    # Preparinng computations of the response matrices
    MMat, tabaMcoef, tabωminωmax = PrepareM(params)

    tabM!(omgval,MMat,tabaMcoef,tabωminωmax,FHT,params)

    if params.VERBOSE>0
        println("LinearResponse.Mode.ComputeModeTables: MMat constructed.")
    end

    # eigenvalue, eigenvector (for basis projection)
    EV, EM = mevXi(MMat)

    return EV, EM
end


"""
    mevXi(tabM)

minimal eigenvalue of M
"""
function mevXi(tabM::AbstractMatrix{ComplexF64})

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

    # output
    return tabeigvals[indEig], tabeigvecs[:,indEig]
end




"""
    GetModeShape(basis,l,Rmin,Rmax,nR,EigenMode)
Function that computes the radial shape of a given mode
-see Eq. (72) of Hamilton+ (2018)



"""
function GetModeShape(basis::AstroBasis.AbstractAstroBasis,
                      Rmin::Float64,Rmax::Float64,
                      nRMode::Int64,
                      EigenMode::Vector,
                      params::LinearParameters)

    # table of R for which the mode is computed
    radiusvals = LinRange(Rmin,Rmax,nRMode)

    # table containing the radial shape of the mode
    eigentype = typeof(EigenMode[1])
    ModeRadius         = zeros(Float64,nRMode)
    ModePotentialShape = zeros(eigentype,nRMode)
    ModeDensityShape   = zeros(eigentype,nRMode)

    (params.VERBOSE > 1) && println("LinearResponse.Mode.GetModeShape: Starting radius loop...")

    # at each radius, compute the shape of the mode response
    for irad=1:nRMode

        # current value of R
        R = radiusvals[irad]

        # initialise the values
        pval = 0.0
        dval = 0.0

        # for each basis element, add the coefficient contribution
        for np=1:params.nradial

            # add the contribution from the basis elements.
            pval += EigenMode[np]*AstroBasis.getUln(basis,params.lharmonic,np-1,R)
            dval += EigenMode[np]*AstroBasis.getDln(basis,params.lharmonic,np-1,R)

        end

        # log contribution in tables
        ModeRadius[irad] = R
        ModePotentialShape[irad] = pval
        ModeDensityShape[irad] = dval

    end

    h5open(ModeFilename(params), "w") do file
        write(file,"ModeRadius",ModeRadius)         # write tabRMode to file
        write(file,"ModePotentialShape",ModePotentialShape) # write tabShapeMode to file
        write(file,"ModeDensityShape",ModeDensityShape) # write tabShapeMode to file
        # Parameters
        WriteParameters(file,params)
    end

    # return just in case we want to do something else
    return ModeRadius,ModePotentialShape,ModeDensityShape
end
