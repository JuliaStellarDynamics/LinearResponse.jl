
"""
    MakeGu(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,n1,n2,np,nq,tabWMat,tabΩ1Ω2Mat,tabAEMat,tabJMat,tabu,Kv,ndim,nradial,ωmin,ωmax,tabvminvmax,lharmonic[,Ω₀])

function to compute G(u)
"""
function MakeGu(ndFdJ::Function,
                n1::Int64,n2::Int64,
                np::Int64,nq::Int64,
                tabWMat::Array{Float64,3},
                tabΩ1Ω2Mat::Array{Float64,3},
                tabELMat::Array{Float64,3},
                tabJMat::Array{Float64,2},
                tabu::Array{Float64},
                Kv::Int64,
                ndim::Int64,
                ωmin::Float64,ωmax::Float64,
                tabvminmax::Array{Float64,2},
                lharmonic::Int64;
                Ω₀::Float64=1.,
                VERBOSE::Int64=0.)

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

    # set up a blank array
    tabGXi  = zeros(Ku)

    for kuval in 1:Ku

        #if (VERBOSE>2)
        #    println("CallAResponse.GFunc.MakeGu: Step $kuval of $Ku.")
        #end

        uval = tabu[kuval]
        vmin, vmax = tabvminmax[kuval,:]

        # determine the step size in v
        deltav = (vmax - vmin)/(Kv)

        res = 0.0 # Initialising the result

        for kvval in 1:Kv
            vval = vmin + deltav*(kvval-0.5)

            ####
            # (u,v) -> (a,e)
            # From WMat computations
            ####
            # (u,v) -> (Ω1,Ω2)
            Ω1,Ω2 = tabΩ1Ω2Mat[kuval,kvval,:]

            # need (E,L): this has some relatively expensive switches
            Eval,Lval = tabELMat[kuval,kvval,1],tabELMat[kuval,kvval,2]

            # compute Jacobians
            #(α,β) -> (u,v).
            # owing to the remapping of Ω, this has an extra 2/(ωmax-ωmin)
            Jacαβ = OrbitalElements.JacαβToUV(n1,n2,ωmin,ωmax,vval)

            #(E,L) -> (α,β): this is the most expensive function here,
            # so we have pre-tabulated it
            JacEL = tabJMat[kuval,kvval]

            #(J) -> (E,L)
            JacJ = (1/Ω1)

            # remove dimensionality from Ω mapping
            dimensionl = (1/Ω₀)

            # get the resonance vector
            ndotΩ = n1*Ω1 + n2*Ω2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotΩ)

            # get tabulated W values for different basis functions np,nq
            Wp = tabWMat[np,kuval,kvval]
            Wq = tabWMat[nq,kuval,kvval]

            if ndim==2
                res += pref*(dimensionl*Jacαβ*JacEL*JacJ*valndFdJ)*Wp*Wq # Local increment in the location (u,v)

            else
                # add in extra Lval from the action-space volume element (Hamilton et al. 2018, eq 30)
                res += pref*Lval*(dimensionl*Jacαβ*JacEL*JacJ*valndFdJ)*Wp*Wq # Local increment in the location (u,v)
            end

        end

        # complete the integration
        res *= deltav
        tabGXi[kuval] = res

    end
    return tabGXi

end


"""
    RunGfunc()

@TO DESCRIBE
"""
function RunGfunc(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                  ndFdJ::Function,
                  wmatdir::String,gfuncdir::String,
                  FHT::FiniteHilbertTransform.FHTtype,
                  Kv::Int64,Kw::Int64,
                  basis::AstroBasis.Basis_type,
                  lharmonic::Int64,
                  n1max::Int64,
                  Ω₀::Float64,
                  modelname::String,dfname::String,
                  rmin::Float64,rmax::Float64;
                  VERBOSE::Int64=0,
                  OVERWRITE::Bool=false)

    # get basis parameters
    ndim, nradial, rb = basis.dimension, basis.nmax, basis.rb

    # Check directory names
    checkdirs = CheckConfigurationDirectories(wmatdir=wmatdir,gfuncdir=gfuncdir)
    if checkdirs < 0
        return 0
    end

    # Integration points
    tabu, Ku = FHT.tabu, FHT.Ku

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    # Frequency cuts associated to [rmin,rmax]
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)

    if VERBOSE >= 0
        println("CallAResponse.GFunc.RunGfunc: Considering $nbResVec resonances.")
    end

    Threads.@threads for i = 1:nbResVec
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        if VERBOSE > 0
            println("CallAResponse.GFunc.RunGfunc: Starting on ($n1,$n2).")
        end

        outputfilename = GFuncFilename(gfuncdir,modelname,dfname,lharmonic,n1,n2,rb,Ku,Kv)
        if isfile(outputfilename)
            if OVERWRITE
                if VERBOSE > 0
                    println("CallAResponse.GFunc.RunGfunc: file already exists for step $i of $nbResVec, ($n1,$n2): recomputing and overwritting.")
                end
            else
                if VERBOSE > 0
                    println("CallAResponse.GFunc.RunGfunc: file already exists for step $i of $nbResVec, ($n1,$n2): no computation.")
                end
                continue
            end
        end

        # load a value of tabWmat, plus (a,e) values
        filename   = WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw)
        file       = h5open(filename,"r")
        Wtab       = read(file,"wmat")
        Ω1Ω2tab    = read(file,"Omgmat")
        AEtab      = read(file,"AEmat")
        Jtab       = read(file,"jELABmat")
        vminmaxtab = read(file,"tabvminmax")
        ωminmax    = read(file,"omgminmax")
        close(file)

        # create the (E,L) table in advance: this does not change for (np,nq)
        ELtab = zeros(Ku,Kv,2)
        for kuval in 1:Ku
            for kvval in 1:Kv
                a,e = AEtab[kuval,kvval,:]

                # need (E,L): this has some relatively expensive switches
                ELtab[kuval,kvval,1],ELtab[kuval,kvval,2] = OrbitalElements.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=0.001)
            end
        end

        # Extremal n.Ω
        ωmin, ωmax = ωminmax[1], ωminmax[2]

        # print the size of the found files if the first processor
        if i==0
            println("CallAResponse.GFunc.RunGfunc: Found nradial=$nradial,Ku=$Ku,Kv=$Kv")
        end

        # need to loop through all combos of np and nq to make the full matrix.
        h5open(outputfilename, "w") do file

            # loop through all basis function combinations
            # can we just do an upper half calculation here?
            for np = 1:nradial
                for nq = 1:nradial

                    if (np==1) & (nq==1)
                        @time tabGXi = MakeGu(ndFdJ,n1,n2,np,nq,
                                              Wtab,Ω1Ω2tab,ELtab,Jtab,
                                              tabu,Kv,ndim,
                                              ωmin,ωmax,
                                              vminmaxtab,
                                              lharmonic,Ω₀=Ω₀,VERBOSE=VERBOSE)
                    else
                        tabGXi = MakeGu(ndFdJ,n1,n2,np,nq,
                                        Wtab,Ω1Ω2tab,ELtab,Jtab,
                                        tabu,Kv,ndim,
                                        ωmin,ωmax,
                                        vminmaxtab,
                                        lharmonic,Ω₀=Ω₀,VERBOSE=VERBOSE)
                    end

                    write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)
                end
            end
        end
    end
end
