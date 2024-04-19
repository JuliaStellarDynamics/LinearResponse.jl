
########################################################################
#
# Compute G(u) = \int dv G(u,v) for one resonance number
#
########################################################################

"""
    MakeGu(ndFdJ,n1,n2,Wdata,tabu[,params])

function to compute G(u)
"""
function MakeGu(distributionfunction::DistributionFunction,
                n1::Int64,n2::Int64,
                Wdata::FourierTransformedBasisData,
                tabu::Vector{Float64},
                params::LinearParameters)

    @assert length(tabu) == params.Ku "Imcompatible tabu length and parameters"

    # calculate the prefactor based on the dimensionality (defaults to 3d)
    if params.dimension==2
        # 2d prefactor, see Fouvry et al. 2015
        # ATTENTION : Landau prescription for G(u) / (u - ω) not G(u) / (ω - u)
        #             Hence the minus sign.
        pref = - (2.0*pi)^(2)
    else
        # 3d prefactor, see Hamilton et al. 2018
        CMatrix = getCMatrix(params.lharmonic)
        pref    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,params.lharmonic,n2)^(2)/(2.0*params.lharmonic+1.0)
    end

    # get basic parameters
    Ku      = length(tabu)
    nradial = size(Wdata.tabW)[1]

    # Number of basis elements
    nradial = size(Wdata.tabW)[1]

    # set up a blank array
    tabGXi  = zeros(Float64,nradial,nradial,Ku)

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
            valndFdJ = ndFdJ((Eval,Lval),(Ω1,Ω2),resonance,distributionfunction)

            # Common part of the integrand (to every np,nq)
            # True volume element (including Jacobians and normalization factors 
            # from change of coordinates (Jr,L) ↦ (u,vp)) : 
            #       (δvp/Ω₀) * (dv/dvp) * 2/(ωmax-ωmin) * |∂(α,β)/∂(u,v)| * |∂(E,L)/∂(α,β)| * |∂(Jr,L)/∂(E,L)|
            # times the prefactor and n⋅∂F/∂J (integrand in (Jr,L)-space)
            integrand = pref * δvol * Jacv * RenormalizedJacαβ * JacEL * JacJ * valndFdJ
            # In 3D, volume element to add
            integrand *= (params.dimension == 3) ? Lval : 1.0

            # Adding step contribution to every element
            for np = 1:nradial
                Wp = Wdata.tabW[np,kvval,kuval]

                # CAUTION: np = nq elements apart (to prevent counting twice)
                tabGXi[np,np,kuval] += integrand * Wp * Wp
                for nq = (np+1):nradial
                    # get tabulated W values for different basis functions np,nq
                    fullin = integrand * Wp * Wdata.tabW[nq,kvval,kuval]

                    # Adding the contribution
                    tabGXi[nq,np,kuval] += fullin
                    tabGXi[np,nq,kuval] += fullin

                end
            end
        end

    end

    return tabGXi
end

########################################################################
#
# Wrapper compute G(u) for all resonance number
#
########################################################################

"""
    RunGfunc(ndFdJ,FHT,basis[,params])

@TO DESCRIBE
"""
function RunGfunc(distributionfunction::DistributionFunction,
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
        tabGXi = MakeGu(distributionfunction,n1,n2,Wdata,FHT.tabu,params)

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
