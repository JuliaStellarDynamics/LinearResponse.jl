
########################################################################
#
# Compute G(u) = \int dv G(u,v) for one resonance number
#
########################################################################


"""
    _ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},res::Resonance, distributionfunction::DistributionFunction)

Distribution function derivative for `distributionfunction` for a given `E`,`L`.

    TODO: make frequencies an optional input
"""
function _ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::EnergyAngularMomentumDF)

        DFDEval,DFDLval = gradient(EL,df)

        # still allow these to be an optional argument
        Ω1,Ω2   = ΩΩ
        n1,n2   = resonance.number[1],resonance.number[2]
        ndotΩ   = n1*Ω1 + n2*Ω2

        return DFDEval * ndotΩ + DFDLval * n2

    end




"""
    MakeGuR(ndFdJ,n1,n2,Wdata,tabu[,params])

function to compute G(u)

    we pick an lmax, and then need to go even l values up to that (no l=zero coupling?)
"""
function MakeGuR(distributionfunction::DistributionFunction,
                n1::Int64,n2::Int64,n3::Int64,
                Wdata::FourierTransformedBasisData,
                tabu::Vector{Float64},
                params::LinearParameters)

    @assert length(tabu) == params.Ku "Imcompatible tabu length and parameters"

    # we've now selected some lmax, called lharmonic
    CMatrix = getCMatrix(params.lharmonic)
    pref0    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,params.lharmonic,n2)^(2)/(2.0*params.lharmonic+1.0)

    # get basic parameters
    Ku      = length(tabu)
    nradial = size(Wdata.tabW)[1]

    # Number of basis elements
    nradial = size(Wdata.tabW)[1]

    # set up a blank array
    tabGXi0  = zeros(Float64,nradial,nradial,Ku)
    tabGXi1a  = zeros(Float64,nradial,nradial,Ku)
    tabGXi1b  = zeros(Float64,nradial,nradial,Ku)

    # determine the step size in vp
    δvp = 1.0/params.Kv

    # remove dimensionality from Ω mapping
    dimensionl = (1.0/Wdata.Ω₀)

    # Integration step volume
    δvol = δvp * dimensionl

    for kuval in 1:params.Ku

        (params.VERBOSE>2) && println("LinearResponse.GFunc.MakeGu: Step $kuval of $(Parmeters.Ku).")

        uval = tabu[kuval]
        vmin, vmax = Wdata.tabvminmax[1,kuval], Wdata.tabvminmax[2,kuval]

        for kvval = 1:params.Kv

            # vp -> v
            vp = δvp*(kvval-0.5)
            vval = v_from_vp(vp, vmin, vmax, n=params.VMAPN)

            ####
            # (u,v) -> (a,e)
            # From WMat computations
            ####
            # (u,v) -> (Ω1,Ω2)
            Ω1,Ω2 = Wdata.tabΩ1Ω2[1,kvval,kuval], Wdata.tabΩ1Ω2[2,kvval,kuval]

            # need (E,L): this has some relatively expensive switches
            Eval,Lval = Wdata.tabEL[1,kvval,kuval], Wdata.tabEL[2,kvval,kuval]
            
            #####
            # compute Jacobians
            #####
            # vp -> v
            _, Jacv = v_from_vp_derivative(vp, vmin, vmax, n=params.VMAPN)
            
            # (u,v) -> (α,β).
            # Renormalized. (2/(ωmax-ωmin) * |∂(α,β)/∂(u,v)|)
            resonance = Resonance(n1,n2,Wdata.ωmin,Wdata.ωmax)
            RenormalizedJacαβ = (2/(Wdata.ωmax-Wdata.ωmin)) * uv_to_αβ_jacobian(uval,vval,resonance)

            # (α,β) -> (E,L): this is the most expensive function here,
            # so we have pre-tabulated it
            JacEL = Wdata.tabJ[kvval,kuval]

            # (E,L) -> (Jr,L)
            JacJ = (1/Ω1)

            # get the resonance vector
            ndotΩ = n1*Ω1 + n2*Ω2
            # compute dF/dJ: call out for value
            #valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotΩ)
            valndFdJ = _ndFdJ((Eval,Lval),(Ω1,Ω2),resonance,distributionfunction)
            valdF = DistributionFunction((Eval,Lval),distributionfunction)

            # Common part of the integrand (to every np,nq)
            # True volume element (including Jacobians and normalization factors 
            # from change of coordinates (Jr,L) ↦ (u,vp)) : 
            #       (δvp/Ω₀) * (dv/dvp) * 2/(ωmax-ωmin) * |∂(α,β)/∂(u,v)| * |∂(E,L)/∂(α,β)| * |∂(Jr,L)/∂(E,L)|
            # times the prefactor and n⋅∂F/∂J (integrand in (Jr,L)-space)

            # compute the jacobian chain: we can probably pull this whole thing out
            jacobianchain = δvol * Jacv * RenormalizedJacαβ * JacEL * JacJ

            # do the non-rotating integrand
            integrand0 = Lval * pref0 * jacobianchain * valndFdJ

            # do the rotating integrand
            integrand1a = Lval * Ifactor(lvalp,lvalq,n2,n3) * jacobianchain * valndFdJ
            integrand1b = Jfactor(lvalp,lvalq,n2,n3) * jacobianchain * valdF


            # Adding step contribution to every element
            for np = 1:nradial
                Wp = Wdata.tabW[np,kvval,kuval]

                # CAUTION: np = nq elements apart (to prevent counting twice)
                tabGXi0[np,np,kuval] += integrand0 * Wp * Wp
                tabGXi1a[np,np,kuval] += integrand1a * Wp * Wp
                tabGXi1b[np,np,kuval] += integrand1b * Wp * Wp
                for nq = (np+1):nradial
                    # get tabulated W values for different basis functions np,nq
                    fullin0 = integrand0 * Wp * Wdata.tabW[nq,kvval,kuval]
                    fullin1a = integrand1a * Wp * Wdata.tabW[nq,kvval,kuval]
                    fullin1b = integrand1b * Wp * Wdata.tabW[nq,kvval,kuval]

                    # Adding the contribution
                    tabGXi0[nq,np,kuval] += fullin0
                    tabGXi0[np,nq,kuval] += fullin0

                    tabGXi1a[nq,np,kuval] += fullin1a
                    tabGXi1a[np,nq,kuval] += fullin1a
                    tabGXi1b[nq,np,kuval] += fullin1b
                    tabGXi1b[np,nq,kuval] += fullin1b

                end
            end
        end

    end

    return tabGXi0,tabXi1a,tabXi1b
end

########################################################################
"""
    RunGfuncR(distributionfunction::DistributionFunction,
              FHT::FiniteHilbertTransform.AbstractFHT,
              params::LinearParameters)

Compute the G(u) function for a given distribution function, finite Hilbert transform, and linear parameters.

# Arguments
- `distributionfunction::DistributionFunction`: The distribution function used to compute G(u).
- `FHT::FiniteHilbertTransform.AbstractFHT`: The finite Hilbert transform used in the computation.
- `params::LinearParameters`: The linear parameters used in the computation.

# Details
This function computes the G(u) function for a given distribution function, finite Hilbert transform, and linear parameters. It iterates over a set of resonances and performs the following steps for each resonance:
1. Checks the directories and FHT values against the parameters.
2. Loads the necessary data from the WMat file.
3. Computes the G(u) values for the resonance.
4. Saves the computed G(u) values in a file.

# Example

"""
function RunGfuncR(distributionfunction::DistributionFunction,
                  FHT::FiniteHilbertTransform.AbstractFHT,
                  params::LinearParameters)

    # check the directories + FHT values against the Parameters
    CheckDirectories(params.wmatdir,params.gfuncdir)
    CheckFHTCompatibility(FHT,params)

    (params.VERBOSE >= 0) && println("LinearResponse.GFunc.RunGfunc: Considering $(params.nbResVec) resonances.")

    Threads.@threads for i = 1:params.nbResVec
        n1, n2 = params.tabResVec[1,i], params.tabResVec[2,i]

        (params.VERBOSE > 0) && println("LinearResponse.GFunc.RunGfunc: Starting on ($n1,$n2).")

        outputfilename = GFuncFilename(n1,n2,params)
        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,params,"LinearResponse.GFunc.RunGfunc: ($n1,$n2) resonance") || continue

        # load a value of tabWmat, plus (a,e) values
        wmatfilename   = WMatFilename(n1,n2,params)
        file       = h5open(wmatfilename,"r")
        # Check if enough basis element in this file (throw error if not)
        (params.nradial <= read(file,"LinearParameters/nradial")) || error("Not enough basis element in WMat file for ($n1,$n2) resonance.")
        # Construct the Wdata structure for the file
        Wdata      = FourierTransformedBasisData(read(file,"omgmin"),read(file,"omgmax"),read(file,"Ω₀"),read(file,"tabvminmax"),                  # Mapping parameters
                                  read(file,"wmat"),                                                                # Basis FT
                                  read(file,"UVmat"),read(file,"Omgmat"),read(file,"AEmat"),read(file,"ELmat"),     # Mappings
                                  read(file,"jELABmat"))                                                            # Jacobians
        close(file)

        # G(u) computation for this resonance number
        tabGXi0,tabXi1a,tabXi1b = MakeGuR(distributionfunction,n1,n2,n3,Wdata,FHT.tabu,params)

        # Saving in file
        h5open(outputfilename, "w") do file
            # Mappings parameters
            write(file, "omgmin",Wdata.ωmin)
            write(file, "omgmax",Wdata.ωmax)
            # G(u)
            write(file,"Gmat",tabGXi)
            # Parameters
            WriteParameters(file,params)
        end
    end
end
