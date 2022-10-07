
"""
    MakeGu(ndFdJ,n1,n2,Wdata,tabu,Kv,ndim,nradial,ωmin,ωmax,tabvminvmax,lharmonic[,Ω₀])

function to compute G(u)
"""
function MakeGu(ndFdJ::Function,
                n1::Int64,n2::Int64,
                Wdata::WMatdata_type,
                tabu::Array{Float64},
                Kv::Int64,
                ndim::Int64,
                lharmonic::Int64;
                Ω₀::Float64=1.,
                VERBOSE::Int64=0)

    # calculate the prefactor based on the dimensionality (defaults to 3d)
    if ndim==2
        # 2d prefactor, see Fouvry et al. 2015
        # ATTENTION : Landau prescription for G(u) / (u - ω) not G(u) / (ω - u)
        #             Hence the minus sign.
        pref = - (2.0*pi)^(2)
    else
        # 3d prefactor, see Hamilton et al. 2018
        CMatrix = getCMatrix(lharmonic)
        pref    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*lharmonic+1.0)
    end

    # get basic parameters
    Ku     = length(tabu)

    # Number of basis elements
    nradial = size(Wdata.tabW)[1]
    # set up a blank array
    tabGXi  = zeros(Float64,nradial,nradial,Ku)
    # get ωmin and ωmax
    ωmin, ωmax = Wdata.ωminmax[1], Wdata.ωminmax[2]

    # determine the step size in vp
    δvp = 1.0/Kv

    # remove dimensionality from Ω mapping
    dimensionl = (1/Ω₀)

    # Integration step volume
    δvol = δvp * dimensionl

    for kuval = 1:Ku

        uval = tabu[kuval]
        vmin, vmax = Wdata.tabvminmax[1,kuval], Wdata.tabvminmax[2,kuval]

        for kvval = 1:Kv
            
            # get the current v value
            vp = δvp*(kvval-0.5)
            vval = g(vp,vmin,vmax,n=2)

            # vp -> v
            Jacvp = dg(vp,vmin,vmax,n=2)

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
            integrand *= (ndim == 3) ? Lval : 1.0

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

"""
    RunGfunc()

@TO DESCRIBE
"""
function RunGfunc(ndFdJ::Function,
                  wmatdir::String,gfuncdir::String,
                  FHT::FiniteHilbertTransform.FHTtype,
                  Kv::Int64,Kw::Int64,
                  basis::AstroBasis.Basis_type,
                  lharmonic::Int64,
                  n1max::Int64,
                  Ω₀::Float64,
                  modelname::String,dfname::String;
                  VERBOSE::Int64=0,
                  OVERWRITE::Bool=false)

    # get basis parameters
    ndim, nradial, rb = basis.dimension, basis.nmax, basis.rb

    # Check directory names
    CheckConfigurationDirectories([wmatdir,gfuncdir]) || (return 0)

    # Integration points
    tabu, Ku = FHT.tabu, FHT.Ku

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    (VERBOSE >= 0) && println("CallAResponse.GFunc.RunGfunc: Considering $nbResVec resonances.")

    Threads.@threads for i = 1:nbResVec
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        (VERBOSE > 1) && println("CallAResponse.GFunc.RunGfunc: Starting on ($n1,$n2).")

        outputfilename = GFuncFilename(gfuncdir,modelname,dfname,lharmonic,n1,n2,rb,Ku,Kv)
        if isfile(outputfilename)
            if OVERWRITE
                (VERBOSE > 0) && println("CallAResponse.GFunc.RunGfunc: file already exists for step $i of $nbResVec, ($n1,$n2): recomputing and overwritting.")
            else
                (VERBOSE > 0) && println("CallAResponse.GFunc.RunGfunc: file already exists for step $i of $nbResVec, ($n1,$n2): no computation.")
                continue
            end
        end

        # load a value of tabWmat, plus (a,e) values
        filename   = WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw)
        file       = h5open(filename,"r")
        Wdata      = WMatdata_type(read(file,"wmat"),                                                               # Basis FT
                                   read(file,"UVmat"),read(file,"Omgmat"),read(file,"AEmat"),read(file,"ELmat"),    # Mappings
                                   read(file,"jELABmat"),                                                           # Jacobians
                                   read(file,"omgminmax"),read(file,"tabvminmax"))                                  # Mapping parameters
        close(file)

        # print the size of the found files if the first processor
        if i==1
            println("CallAResponse.GFunc.RunGfunc: Found nradial=$nradial,Ku=$Ku,Kv=$Kv")
        end

        # G(u) computation for this resonance number
        tabGXi = MakeGu(ndFdJ,n1,n2,
                        Wdata,
                        tabu,Kv,ndim,
                        lharmonic,Ω₀=Ω₀,VERBOSE=VERBOSE)
        # Saving in file
        h5open(outputfilename, "w") do file
        write(file,"Gmat",tabGXi)
        end
    end
end
