



"""
for a single omega, compute the shape of the mode

"""
function ComputeModeTables(inputfile,
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
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # get all Legendre weights
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(K_u)

    # make the (np,nq) vectors that we need to evaluate
    tab_npnq = makeTabnpnq(nradial)

    # make the decomposition coefficients a_k
    MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE)

    # load aXi values
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial)
    println("CallAResponse.Xi.FindZeroCrossing: tabaMcoef loaded.")

    # struct for D_k(omega) computation
    struct_tabLeglist = FiniteHilbertTransform.struct_tabLeg_create(K_u)

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1
    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    tabM!(omgval,MMat,tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist,dψ,d2ψ,nradial,Omega0)
    println("CallAResponse.Mode.ComputeModeTables: MMat constructed.")

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
Function that computes the radial shape of a given mode
-see Eq. (72) of Hamilton+ (2018)
"""
function GetModeShape(inputfile::String,
                      Rmin::Float64,Rmax::Float64,
                      nRMode::Int64,
                      tabEigenMode::Array{Float64})

    # load model parameters: this exposes basis
    include(inputfile)

    # prep the basis
    AstroBasis.fill_prefactors!(basis)

    # get the step distance of the array in RMode
    deltaRMode = (Rmax - Rmin)/(nRMode - 1)

    # table of R for which the mode is computed
    tabRMode = collect(Rmin:deltaRMode:Rmax)

    # table containing the radial shape of the mode
    tabShapeMode = zeros(Float64,nRMode)

    println("CallAResponse.Mode.GetModeShape: Starting radius loop...")

    # at each radius, compute the shape of the mode response
    for irad=1:nRMode

        # current value of R
        R = tabRMode[irad]

        # initialise the value
        val = 0.0

        # for each basis element, add the coefficient contribution
        for np=1:nradial

            # add the contribution from the basis elements.
            val += tabEigenMode[np]*AstroBasis.getUln(basis,lharmonic,np-1,R)

        end

        # log contribution in table
        tabShapeMode[irad] = val
    end

    h5open(mode_filename(modedir,modelname,lharmonic,n1max,K_u), "w") do file
        write(file,"tabRMode",tabRMode)         # write tabRMode to file
        write(file,"tabShapeMode",tabShapeMode) # write tabShapeMode to file
    end

    # return just in case we want to do something else
    return tabRMode,tabShapeMode

end
