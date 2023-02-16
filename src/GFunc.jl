
########################################################################
#
# Compute G(u) = \int dv G(u,v) for one resonance number
#
########################################################################

"""
    MakeGu(ndFdJ,n1,n2,Wdata,tabu,params)

function to compute G(u)
"""
function MakeGu(ndFdJ::Function,
                n1::Int64,n2::Int64,
                Wdata::WMatdataType,
                tabu::Array{Float64},
                Parameters::ResponseParameters)

    # calculate the prefactor based on the dimensionality (defaults to 3d)
    if Parameters.ndim==2
        # 2d prefactor, see Fouvry et al. 2015
        # ATTENTION : Landau prescription for G(u) / (u - ω) not G(u) / (ω - u)
        #             Hence the minus sign.
        pref = - (2.0*pi)^(2)
    else
        # 3d prefactor, see Hamilton et al. 2018
        CMatrix = getCMatrix(Parameters.lharmonic)
        pref    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,Parameters.lharmonic,n2)^(2)/(2.0*Parameters.lharmonic+1.0)
    end

    # get basic parameters
    Ku      = length(tabu)
    nradial = size(Wdata.tabW)[1]

    # Number of basis elements
    nradial = size(Wdata.tabW)[1]

    # set up a blank array
    tabGXi  = zeros(Float64,nradial,nradial,Ku)

    # get ωmin and ωmax
    ωmin, ωmax = Wdata.ωmin, Wdata.ωmax

    # determine the step size in vp
    δvp = 1.0/Parameters.Kv

    # remove dimensionality from Ω mapping
    dimensionl = (1.0/Parameters.OEparams.Ω₀)

    # Integration step volume
    δvol = δvp * dimensionl

    for kuval in 1:Parameters.Ku

        (Parameters.VERBOSE>2) && println("CallAResponse.GFunc.MakeGu: Step $kuval of $(Parmeters.Ku).")

        uval = tabu[kuval]
        vmin, vmax = Wdata.tabvminmax[1,kuval], Wdata.tabvminmax[2,kuval]

        for kvval = 1:Parameters.Kv

            # get the current v value
            vp = δvp*(kvval-0.5)
            vval = vFromvp(vp,vmin,vmax,Parameters.VMAPN)

            # vp -> v
            Jacvp = DvDvp(vp,vmin,vmax,Parameters.VMAPN)

            ####
            # (u,v) -> (a,e)
            # From WMat computations
            ####
            # (u,v) -> (Ω1,Ω2)
            Ω1,Ω2 = Wdata.tabΩ1Ω2[1,kvval,kuval], Wdata.tabΩ1Ω2[2,kvval,kuval]

            # need (E,L): this has some relatively expensive switches
            Eval,Lval = Wdata.tabEL[1,kvval,kuval], Wdata.tabEL[2,kvval,kuval]

            # compute Jacobians
            # (α,β) -> (u,v).
            # owing to the remapping of Ω, this has an extra 2/(ωmax-ωmin)
            Jacαβ = OrbitalElements.JacαβToUV(n1,n2,ωmin,ωmax,vval)

            # (E,L) -> (α,β): this is the most expensive function here,
            # so we have pre-tabulated it
            JacEL = Wdata.tabJ[kvval,kuval]

            #(J) -> (E,L)
            JacJ = (1/Ω1)


            # get the resonance vector
            ndotΩ = n1*Ω1 + n2*Ω2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotΩ)


            # Common part of the integrand (to every np,nq)
            integrand = pref*δvol*Jacvp*Jacαβ*JacEL*JacJ*valndFdJ
            # In 3D, volume element to add
            integrand *= (Parameters.ndim == 3) ? Lval : 1.0

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
    RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,FHT,basis,Parameters)

@TO DESCRIBE
"""
function RunGfunc(ndFdJ::Function,
                  FHT::FiniteHilbertTransform.FHTtype,
                  Parameters::ResponseParameters)

    # Check directory names
    CheckConfigurationDirectories([Parameters.wmatdir,Parameters.gfuncdir]) || (return 0)

    (Parameters.VERBOSE >= 0) && println("CallAResponse.GFunc.RunGfunc: Considering $(Parameters.nbResVec) resonances.")

    Threads.@threads for i = 1:Parameters.nbResVec
        n1, n2 = Parameters.tabResVec[1,i], Parameters.tabResVec[2,i]

        (Parameters.VERBOSE > 0) && println("CallAResponse.GFunc.RunGfunc: Starting on ($n1,$n2).")

        outputfilename = GFuncFilename(n1,n2,Parameters)
        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,Parameters,"CallAResponse.GFunc.RunGfunc: ($n1,$n2) resonance") || continue

        # load a value of tabWmat, plus (a,e) values
        wmatfilename   = WMatFilename(n1,n2,Parameters)
        file       = h5open(wmatfilename,"r")
        # Check if enough basis element in this file (throw error if not)
        (Parameters.nradial <= read(file,"ResponseParameters/nradial")) || error("Not enough basis element in WMat file for ($n1,$n2) resonance.")
        # Construct the Wdata structure for the file
        Wdata      = WMatdataType(read(file,"omgmin"),read(file,"omgmax"),read(file,"tabvminmax"),                  # Mapping parameters
                                  read(file,"wmat"),                                                                # Basis FT
                                  read(file,"UVmat"),read(file,"Omgmat"),read(file,"AEmat"),read(file,"ELmat"),     # Mappings
                                  read(file,"jELABmat"))                                                            # Jacobians
        close(file)

        # G(u) computation for this resonance number
        tabGXi = MakeGu(ndFdJ,n1,n2,Wdata,FHT.tabu,Parameters)

        # Saving in file
        h5open(outputfilename, "w") do file
            # Mappings parameters
            write(file, "omgmin",Wdata.ωmin)
            write(file, "omgmax",Wdata.ωmax)
            # G(u)
            write(file,"Gmat",tabGXi)
            # Parameters
            WriteParameters(file,Parameters)
        end
    end
end
