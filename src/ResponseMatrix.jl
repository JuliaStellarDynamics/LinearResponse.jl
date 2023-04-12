"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

########################################################
# Compute response matrix at a given (complex) frequency
########################################################

"""
    tabM!(ω,tabM,tabaMcoef,tabωminωmax,FHT,params)

computes the response matrix M(ω) for a given COMPLEX frequency ω in physical units, i.e. not (yet) rescaled by 1/Ω0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabM!(ω::ComplexF64,
               tabM::AbstractMatrix{ComplexF64},
               tabaMcoef::Array{Float64,4},
               tabωminωmax::Matrix{Float64},
               FHT::FiniteHilbertTransform.AbstractFHT,
               params::LinearParameters)

    # get dimensions from the relevant tables
    nbResVec, tabResVec = params.nbResVec, params.tabResVec
    KuTruncation = params.KuTruncation
    nradial  = params.nradial
    VERBOSE  = params.VERBOSE
    Ku       = FHT.Ku

    # Assertions before computation
    CheckFHTCompatibility(FHT,params)
    CheckMShape(tabM;nradial=nradial) 
    CheckaMcoefShape(tabaMcoef;nradial=nradial,Ku=Ku,nbResVec=nbResVec) 

    if KuTruncation < Ku
        Ku = KuTruncation
        (VERBOSE > 2) && println("LinearResponse.Xi.tabM!: truncating Ku series from $(FHT.Ku) to $Ku.")
    end

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    # Rescale to get dimensionless frequency
    ωnodim = ω/params.Orbitalparams.Ω₀

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nres = 1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        # get ωmin and ωmax values
        ωmin, ωmax = tabωminωmax[1,nres], tabωminωmax[2,nres]

        # get the rescaled frequency
        ϖ = OrbitalElements.Getϖ(ωnodim,ωmin,ωmax)

        # get the integration values
        FiniteHilbertTransform.GettabD!(ϖ,FHT)

        # mame of the array where the D_k(w) are stored
        tabD = FHT.tabDLeg

        # loop over the basis indices to consider
        for np = 1:nradial
            for nq = np:nradial

                res = 0.0 + 0.0*im

                # loop over the Legendre functions to add all contributions
                for k=1:Ku

                    # hard check for nans
                    @inbounds val = tabaMcoef[k,nq,np,nres]*tabD[k]
                    if !isnan(val)
                        res += val
                    else
                        (k==1) && (VERBOSE>1) && println("LinearResponse.Xi.tabM!: NaN found for n=($n1,$n2), npnq=($np,$nq), k=$k")
                    end


                end

                # fill the full M matrix:
                # as tab_npnq is the upper triangular matrix (with the diagonal),
                # we need to duplicate for symmetries

                # fill in the element (np,nq)
                @inbounds tabM[np,nq] += res

                # fill in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.
                @inbounds tabM[nq,np] += res
            end
        end # basis index loop

    end # resonance loop

    # contributions were added twice for diagonal elements: reset
    for np=1:nradial
        @inbounds tabM[np,np] *= 0.5
    end
end


########################################
# To prepare response matrix computation
########################################

"""

prepare computations for Linear Response (Determinant, eigenvalues, mode, full matrices computations...)
for multiple threads.
"""
function PrepareM(n::Int64,
                  FHT::FiniteHilbertTransform.AbstractFHT,
                  params::LinearParameters)

    
    @assert n > 0

    # Check directory names
    CheckDirectories(params.axidir,params.modedir)

    # allocate memory for FHT structs
    FHTlist = [deepcopy(FHT) for k=1:n]

    # allocate memory for the response matrices M
    nradial = params.nradial
    tabMlist = [zeros(ComplexF64,nradial,nradial) for k=1:n]

    # load aXi values
    tabaMcoef, tabωminωmax = StageAXi(params)

    (params.VERBOSE >= 1) && println("LinearResponse.Xi.RunM: tabaMcoef loaded.")

    return tabMlist, tabaMcoef, tabωminωmax, FHTlist
end

"""

prepare computations for Linear Response (Determinant, eigenvalues, mode, full matrices computations...)
for single thread.
"""
function PrepareM(params::LinearParameters)

    # check the directories
    CheckDirectories(params.axidir,params.modedir)

    # allocate memory for the response matrice M
    nradial = params.nradial
    MMat = zeros(ComplexF64,nradial,nradial)

    # load aXi values
    tabaMcoef, tabωminωmax = StageAXi(params)

    (params.VERBOSE >= 1) && println("LinearResponse.Xi.RunM: tabaMcoef loaded.")

    return MMat, tabaMcoef, tabωminωmax
end