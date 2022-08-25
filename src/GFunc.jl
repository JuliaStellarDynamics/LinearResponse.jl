
"""
function to compute G(u)

@ATTENTION, the dimensionality (e.g. 2d vs 3d) is now encoded in 'ndim'.

"""
function MakeGu(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                ndFdJ::Function,
                n1::Int64,n2::Int64,
                np::Int64,nq::Int64,
                tabWMat::Array{Float64},
                tabaMat::Array{Float64},
                tabeMat::Array{Float64},
                tabJMat::Array{Float64},
                Kuvals::Matrix{Float64},
                K_v::Int64,nradial::Int64,
                ωmin::Float64,ωmax::Float64,
                vminarr::Array{Float64},vmaxarr::Array{Float64},
                lharmonic::Int64;
                ndim::Int64,
                Ω0::Float64=1.)

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
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    for kuval in 1:K_u

        uval = Kuvals[kuval]
        vmin = vminarr[kuval]
        vmax = vmaxarr[kuval]

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        res = 0.0 # Initialising the result

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            α,β = OrbitalElements.AlphaBetaFromUV(uval,vval,n1,n2,ωmin,ωmax)

            Ω1,Ω2 = α*Ω0,α*β*Ω0

            # convert from Ω1,Ω2 to (a,e): using a tabled value
            a,e = tabaMat[kuval,kvval],tabeMat[kuval,kvval]

            # need (E,L): this has some relatively expensive switches
            Eval,Lval = OrbitalElements.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=0.001)

            # compute Jacobians
            #(α,β) -> (u,v).
            # owing to the remapping of Ω, this has an extra 2/(ωmax-ωmin)
            Jacαβ = OrbitalElements.JacAlphaBetaToUV(n1,n2,ωmin,ωmax,vval)

            #(E,L) -> (α,β): this is the most expensive function here
            #JacEL        = OrbitalElements.JacEL_to_alphabeta(alpha,beta)
            #JacEL        = OrbitalElements.JacELToAlphaBetaAE(a,e,ψ,dψ,d2ψ,Ω0)
            JacEL = tabJMat[kuval,kvval]
            #JacEL = OrbitalElements.JacELToAlphaBetaAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,Ω0=Ω0)

            #(J) -> (E,L)
            JacJ = (1/Ω1)

            # remove dimensionality from Ω mapping
            dimensionl = (1/Ω0)


            # get the resonance vector
            ndotΩ = n1*Ω1 + n2*Ω2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotΩ)

            # get tabulated W values for different basis functions np,nq
            Wp = tabWMat[np,kuval,kvval]
            Wq = tabWMat[nq,kuval,kvval]

            # todo: make this block @static
            #=
            # do a nan check?
            nancheck = false
            if (nancheck)
                tmp = pref*Lval*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq

                if isnan(tmp)
                    println(Jacalphabeta," ",JacEL," ",pref," ",(Lval/Ω1)," ",valndFdJ," ",Wp," ",Wq)
                end
            end
            =#

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
    RunGfunc(inputfile)

"""
function RunGfunc(inputfile::String)

    # bring in the defined parameters
    LoadConfiguration(inputfile)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(wmatdir=wmatdir,gfuncdir=gfuncdir)
    if checkdirs < 0
        return 0
    end

    # prep for Legendre integration
    tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
    tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

    # make a function for the circular frequency relationship:
    #  only needs to happen once, redefinition might kill us
    #  could this move completely outside the function to help the compiler?
    βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω0,rmax=1.0e6)

    # compute the number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    println("CallAResponse.GFuncIsochrone.RunGfunc: Considering $nbResVec resonances.")

    Threads.@threads for i = 1:nbResVec
        n1,n2 = tabResVec[1,i],tabResVec[2,i]
        println("CallAResponse.GFuncIsochrone.RunGfunc: Starting on ($n1,$n2).")

        # load a value of tabWmat, plus (a,e) values
        filename = wmat_filename(wmatdir,modelname,lharmonic,n1,n2,rb)
        file = h5open(filename,"r")
        Wtab = read(file,"wmat")
        atab = read(file,"amat")
        etab = read(file,"emat")
        Jtab = read(file,"jELABmat")
        nradial,K_u,K_v = size(Wtab)

        # print the size of the found files if the first processor
        if i==0
            println("CallAResponse.GFuncIsochrone.RunGfunc: Found nradial=$nradial,K_u=$K_u,K_v=$K_v")
        end

        outputfilename = GFuncFilename(gfuncdir,modelname,dfname,lharmonic,n1,n2,K_u,rb)
        if isfile(outputfilename)
            println("CallAResponse.GFuncIsochrone.RunGfunc: file already exists for step $i of $nbResVec, ($n1,$n2).")
            continue
        end

        # compute the frequency scaling factors for this resonance
        ωmin,ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,10000.,Ω0)

        vbound = OrbitalElements.FindVbound(n1,n2,dψ,d2ψ,10000.,Ω0)

        # for some threading reason, make sure K_u is defined here
        K_u = length(tabwGLquad)

        # loop through once and design a v array for min, max
        vminarr,vmaxarr = zeros(K_u),zeros(K_u)
        for uval = 1:K_u
           vminarr[uval],vmaxarr[uval] = OrbitalElements.FindVminVmax(tabuGLquad[uval],ωmin,ωmax,n1,n2,vbound,βc)
        end

        # need to loop through all combos of np and nq to make the full matrix.
        h5open(outputfilename, "w") do file

            # loop through all basis function combinations
            # can we just do an upper half calculation here?
            for np = 1:nradial
                for nq = 1:nradial

                    if (np==1) & (nq==1)
                        @time tabGXi = MakeGu(ψ,dψ,d2ψ,d3ψ,d4ψ,
                                              ndFdJ,n1,n2,np,nq,
                                              Wtab,atab,etab,Jtab,
                                              tabuGLquad,K_v,nradial,
                                              ωmin,ωmax,
                                              vminarr,vmaxarr,
                                              lharmonic,ndim=ndim,Ω0=Ω0)
                    else
                        tabGXi = MakeGu(ψ,dψ,d2ψ,d3ψ,d4ψ,
                                        ndFdJ,n1,n2,np,nq,
                                        Wtab,atab,etab,Jtab,
                                        tabuGLquad,K_v,nradial,
                                        ωmin,ωmax,
                                        vminarr,vmaxarr,
                                        lharmonic,ndim=ndim,Ω0=Ω0)
                    end

                    #

                    #=
                    sumG = sum(tabGXi)
                    if (np>-100) & (nq>-100)
                        if isnan(sumG)
                            println("Gfunc.jl: NaN for n1=$n1, n2=$n2.")
                        else
                            println("GFunc.jl: np=$np, nq=$nq, sumG=$sumG.")
                        end
                    end
                    =#
                    write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)
                end
            end
        end
    end

end
