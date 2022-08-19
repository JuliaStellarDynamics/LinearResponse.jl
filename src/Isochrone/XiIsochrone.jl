"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

using LinearAlgebra


"""tabM!(omg,tabM,tabaMcoef,tabResVec,tabnpnq,struct_tabLeg,LINEAR,Omega0)
Function that computes the response matrix Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Omega0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabMIsochrone!(omg::Complex{Float64},
                        tabM::Array{Complex{Float64},2},
                        tabaMcoef::Array{Float64,4},
                        tabResVec::Matrix{Int64},
                        tab_npnq::Matrix{Int64},
                        struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
                        nradial::Int64,
                        LINEAR::String="unstable",
                        Omega0::Float64=1.0)

    # get dimensions from the relevant tables
    nb_npnq  = size(tab_npnq)[2]
    nbResVec = size(tabResVec)[2]
    K_u      = size(tabaMcoef)[4]

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    tabMtmp = zeros(Complex{Float64},nradial,nradial)

    #println("CallAResponse.Xi.tabM!: loop dimensions npnq=$nb_npnq, nResVec=$nbResVec, K_u=$K_u.")

    # set up to output the tabM values
    #h5open(MFilename(modedir,modelname,dfname,lharmonic,K_u), "w") do file

        # loop over the resonances: no threading here because we parallelise over frequencies
        for nResVec=1:nbResVec

            # get current resonance numbers (n1,n2)
            n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec]

            # Rescale to get dimensionless frequency
            omg_nodim = omg/Omega0

            # get the rescaled frequency
            varpi = OrbitalElements.GetVarpiIsochrone(omg_nodim,n1,n2)

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
                    res += tabaMcoef[nResVec,np,nq,k]*tabDLeg[k]
                end

                # fill the full M matrix:
                # as tab_npnq is the upper triangular matrix (with the diagonal),
                # we need to duplicate for symmetries

                # fill in the element (np,nq)
                tabM[np,nq] += res

                # fill in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.
                tabM[nq,np] += res

                tabMtmp[np,nq] = res
                tabMtmp[nq,np] = res

            end # basis index loop

            # do the actual writing
            #write(file, "tabMn1"*string(n1)*"n2"*string(n2),tabMtmp)

        end # resonance loop

        # contributions were added twice for diagonal elements: reset
        for np=1:nradial
            tabM[np,np] *= 0.5
        end

    #end # the M output file

end




function RunMIsochrone(inputfile::String,
                       omglist::Array{Complex{Float64}};
                       VERBOSE::Int64=0)

    # need some sort of 'if' for whether this already exists
    include(inputfile)

    nomglist = length(omglist)

    #####
    # Check directories names
    #####
    if !(isdir(gfuncdir) && isdir(modedir))
        error("CallAResponse.Xi.RunMIsochrone: gfuncdir or modedir not found.")
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
    MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,OVERWRITE=false)

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
    tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial)
    println("CallAResponse.Xi.RunMIsochrone: tabaMcoef loaded.")

    println("CallAResponse.Xi.RunMIsochrone: Starting frequency analysis, using $LINEAR integration.")
    println("CallAResponse.Xi.RunMIsochrone: computing $nomglist frequency values.")

    # loop through all frequencies
    Threads.@threads for i = 1:nomg

        k = Threads.threadid()

        if i==1
            @time tabMIsochrone!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist[k],nradial,LINEAR,Omega0)
        else
            tabMIsochrone!(omglist[i],tabMlist[k],tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist[k],nradial,LINEAR,Omega0)
        end

        tabdetXi[i] = detXi(IMatlist[k],tabMlist[k])

    end

    WriteDeterminant(det_filename(modedir,modelname,dfname,lharmonic,n1max,K_u),omglist,tabdetXi)

    return tabdetXi
end
