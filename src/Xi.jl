"""

@IMPROVE: what is the best way to pass the format for reading in Gfunctions? (might want to consider a new Gfunc format)

"""

using LinearAlgebra


"""makeaMCoefficients!(tabaMcoef,tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,basedir)

function to make the decomposition coefficients "a" of the response matrix M
"""
function makeaMCoefficients!(tabaMcoef::Array{Float64,4},
                             tabResVec::Matrix{Int64},
                             tabnpnq::Matrix{Int64},
                             tabwGLquad::Vector{Float64},
                             tabPGLquad::Matrix{Float64},
                             tabINVcGLquad::Vector{Float64},
                             gfuncdir::String,
                             modelname::String,
                             lharmonic::Int64)

    # get relevant sizes
    K_u      = size(tabaMcoef)[4]
    nb_npnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    # loop through all resonance vectors
    for nresvec in 1:nbResVec
        n1,n2 = tabResVec[1,nresvec],tabResVec[2,nresvec]

        # open the resonance file
        filename = gfunc_filename(gfuncdir,modelname,lharmonic,n1,n2,K_u)
        file = h5open(filename,"r")

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)

            # open the correct resonance vector
            # read in the correct G(u) function
            tabGXi = read(file,"GXinp"*string(np)*"nq"*string(nq))

            # DEBUG nan values in tabGXi
            #sumG = sum(tabGXi)
            #println(println("n1=$n1, n2=$n2, np=$np, nq=$nq, sumG=$sumG."))

            # Loop over the Legendre functions
            for k=1:K_u

                res = 0.0 # Initialisation of the result

                for i=1:K_u # Loop over the G-L nodes
                    w = tabwGLquad[i] # Current weight
                    G = tabGXi[i] # Current value of G[u_i]
                    P = tabPGLquad[k,i] # Current value of P_k. ATTENTION, to the order of the arguments.
                    res += w*G*P # Update of the sum
                end

                res *= tabINVcGLquad[k] # Multiplying by the Legendre prefactor.

                # populate the symmetric matrix
                tabaMcoef[nresvec,np,nq,k] = res # Element (np,nq)
                tabaMcoef[nresvec,nq,np,k] = res # Element (nq,np). If np=nq overwrite.

            end

            #h5write("xifunc/tabaXI_n1_"*string(n1)*"_n2_"*string(n2)*"_np_"*string(np)*"_nq_"*string(nq)*"."*string(K_u)*".h5","tabGXi",tabaXi[nresvec][np,nq])
        end
    end
end


"""tabM!(omg,tabM,tabaMcoef,tabResVec,tabnpnq,struct_tabLeg,dpotential,ddpotential,LINEAR,Omega0)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Omega0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabM!(omg::Complex{Float64},
                tabM::Array{Complex{Float64},2},
                tabaMcoef::Array{Float64,4},
                tabResVec::Matrix{Int64},
                tab_npnq::Matrix{Int64},
                struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
                dpotential::Function,
                ddpotential::Function,
                nradial::Int64,
                LINEAR::String="unstable",
                Omega0::Float64=1.0)

    # get dimensions from the relevant tables
    nb_npnq  = size(tab_npnq)[2]
    nbResVec = size(tabResVec)[2]
    K_u      = size(tabaMcoef)[4]

    fill!(tabM,0.0 + 0.0*im) # Initialising the array to 0.
    #####
    for nResVec=1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec] # Current resonance (n1,n2)

        # Rescale to get the dimensionless frequency
        omg_nodim = omg/Omega0
        varpi = OrbitalElements.get_varpi(omg_nodim,n1,n2,dpotential,ddpotential,Ω₀=Omega0) # Getting the rescaled frequency

        # get the Legendre integration values
        PerturbPlasma.get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR)
        #PerturbPlasma.tabLeg!_UNSTABLE(varpi,K_u,struct_tabLeg)
        tabDLeg = struct_tabLeg.tabDLeg # Name of the array where the D_k(w) are stored

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tab_npnq[1,i_npnq], tab_npnq[2,i_npnq] # Current value of (np,nq)

            # Loop over the Legendre functions to add all contributions
            res = 0.0 + 0.0*im
            for k=1:K_u
                res += tabaMcoef[nResVec,np,nq,k]*tabDLeg[k]
            end

            # fill the full matrix
            tabM[np,nq] += res # Filling in the element (np,nq)
            tabM[nq,np] += res # Filling in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.

        end # basis index loop
        #h5write("xifunc/tabXI_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5","tabXi",tabXi)

    end # resonance loop

    for np=1:nradial
        tabM[np,np] *= 0.5 # Contributions added twice for diagonal elements
    end

end


# make an identity matrix
function makeIMat(nradial::Int64)
    IMat = zeros(Complex{Float64},nradial,nradial) # Static container for the identity matrix
    ##########
    for np=1:nradial # Loop over the radial elements to fill in the identity matrix
        IMat[np,np] = 1.0 + 0.0im # Creating the identity matrix
    end
    return IMat
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
    for np=1:nEig1 # Loop over the number of basis elements
        tabEigenMode[np] = real(tabeigvecs[np,indEig]) # Extracting the eigenvector !! ATTENTION, tabeigvecs is the matrix whose columns are eigenvectors !! ATTENTION, we put a real part
    end

    return tabeigvals[indEig],tabeigvecs[indEig],tabEigenMode # Output

end


"""
    detXi(.......)

determinant of the susceptibility matrix I - M for unknown M.
"""
function detXi(omg::Complex{Float64},
               tabM::Array{Complex{Float64},2},
               tabaMcoef::Array{Float64,4},
               tabResVec::Matrix{Int64},
               tab_npnq::Matrix{Int64},
               struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
               dpotential::Function,
               ddpotential::Function,
               nradial::Int64,
               LINEAR::String="unstable",
               Omega0::Float64=1.0)
    ####
    IMat = makeIMat(nradial)
    tabM!(omg,tabM,tabaMcoef,tabResVec,tab_npnq,struct_tabLeg,dpotential,ddpotential,nradial,LINEAR,Omega0)
    ####
    detXi(IMat,tabM)
end


function RunM(inputfile,
              omglist::Array{Complex{Float64}})

    include(inputfile)

    nomglist = length(omglist)
    println("Xi.jl: computing $nomglist frequency values")

    #####
    # Check directories names
    #####
    if !(isdir(wmatdir) && isdir(gfuncdir))
        error(" wmatdir or gfuncdir not found ")
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
    struct_tabLeglist = [PerturbPlasma.struct_tabLeg_create(K_u) for k=1:Threads.nthreads()]
    # memory for the response matrices M and identity matrices
    tabMlist = [zeros(Complex{Float64},nradial,nradial) for k=1:Threads.nthreads()]
    IMat = makeIMat(nradial)
    IMatlist = [deepcopy(IMat) for k=1:Threads.nthreads()]

    # Containers for determinant and min eigenvalue
    nomg = length(omglist)
    tabdetXi = zeros(Complex{Float64},nomg) # Real part of the determinant at each frequency
    # tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency


    Threads.@threads for i = 1:nomg

        k = Threads.threadid()

        tabM!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist[k],dpotential,ddpotential,nradial,LINEAR,Omega0)

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])
        #tabmevXi[i] = mevXi(IMatlist[k],tabMlist[k])

    end

    return tabdetXi#, tabmevXi
end
