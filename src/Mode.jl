



"""
for a single omega, compute the shape of the mode

"""
function RunShape(inputfile,
                  omgval::Complex{Float64})

    include(inputfile)

    #####
    # Check directories names
    #####
    if !(isdir(gfuncdir) && isdir(modedir))
        error("Mode.jl: gfuncdir or modedir not found ")
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
    tabaMcoef = zeros(Float64,nbResVec,nradial,nradial,K_u)
    makeaMCoefficients!(tabaMcoef,tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,lharmonic)

    # Structs for D_k(omega) computation
    struct_tabLeglist = PerturbPlasma.struct_tabLeg_create(K_u)

    # memory for the response matrices M and identity matrices
    MMat = zeros(Complex{Float64},nradial,nradial)

    # Containers for determinant and min eigenvalue
    nomg = 1#length(omglist)
    tabdetXi = zeros(Float64,nomg) # Real part of the determinant at each frequency
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    tabM!(omgval,MMat,tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist,dpotential,ddpotential,nradial,LINEAR,Omega0)

    # eigenvalue, eigenfunction (eigenvector), eigenmode (for basis projection)
    EV,EF,EM = mevXi(MMat)

    return EV,EF,EM
end


"""
Function that computes the radial shape of a given mode
-see Eq. (72) of Hamilton+ (2018)
"""
function tabShapeMode(inputfile::String,
                      Rmin::Float64,Rmax::Float64,
                      nRMode::Int64,
                      tabEigenMode::Array{Float64})#,
                      #basis::AstroBasis.Basis_type)

    # load model parameters
    include(inputfile)

    # prep the basis
    AstroBasis.fill_prefactors!(basis)

    # get the step distance of the array in RMode
    deltaRMode = (Rmax - Rmin)/(nRMode - 1)

    # table of R for which the mode is computed
    tabRMode = collect(Rmin:deltaRMode:Rmax)

    # table containing the radial shape of the mode
    tabShapeMode = zeros(Float64,nRMode)

    # at each radius, compute the shape of the mode response
    for irad=1:nRMode

        # current value of R
        R = tabRMode[irad]

        # initialise the value
        val = 0.0

        # for each basis element, add the coefficient contribution
        for np=1:nradial-1

            # add the contribution from the basis elements.
            val += tabEigenMode[np]*AstroBasis.getUln(basis,lharmonic,np,R)

        end

        tabShapeMode[irad] = val
    end

    h5open(mode_filename(modedir,modelname,lharmonic,n1max,K_u), "w") do file
        write(file,"tabRMode",tabRMode)         # write tabRMode to file
        write(file,"tabShapeMode",tabShapeMode) # write tabShapeMode to file
    end
    # disable return
    return tabRMode,tabShapeMode

end
