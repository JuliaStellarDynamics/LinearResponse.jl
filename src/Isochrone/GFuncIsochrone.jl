


"""MakeGuIsochrone
function to compute G(u), isochrone specific.

"""
function MakeGuIsochrone(ndFdJ::Function,
                         n1::Int64,n2::Int64,
                         np::Int64,nq::Int64,
                         tabWMat::Array{Float64},
                         tabaMat::Array{Float64},
                         tabeMat::Array{Float64},
                         Kuvals::Matrix{Float64},
                         K_v::Int64,nradial::Int64,
                         ωmin::Float64,ωmax::Float64,
                         vminarr::Array{Float64},vmaxarr::Array{Float64},
                         lharmonic::Int64;
                         ndim::Int64,
                         Omega0::Float64=1.,bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    # calculate the prefactor based on the dimensionality (defaults to 3d)
    if ndim==2
        pref = (2pi)^2
    else
        CMatrix = getCMatrix(lharmonic)
        pref    = -2.0*(2pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2lharmonic+1)
    end

    # get basic parameters
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    for kuval in 1:K_u

        # retrieve integration values
        uval = Kuvals[kuval]
        vmin = vminarr[kuval]
        vmax = vmaxarr[kuval]

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        # initialise the result
        res = 0.0

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.AlphaBetaFromUV(uval,vval,n1,n2,ωmin,ωmax)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e) using isochrone exact version
            sma,ecc = OrbitalElements.IsochroneAEFromOmega1Omega2(omega1,omega2,bc,M,G)

            # get (rp,ra)
            rp,ra = OrbitalElements.RpRafromAE(sma,ecc)

            # need (E,L), use isochrone exact version
            Eval,Lval = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)

            # compute Jacobians
            #(alpha,beta) -> (u,v)
            Jacalphabeta = OrbitalElements.JacalphabetaToUV(n1,n2,ωmin,ωmax,vval)

            #(E,L) -> (alpha,beta): Isochrone analytic
            JacEL        = OrbitalElements.IsochroneJacELtoAlphaBeta(alpha,beta,bc,M,G)

            #(J) -> (E,L)
            JacJ         = (1/omega1)

            # remove dimensionality
            dimensionl   = (1/Omega0)


            # get the resonance vector
            ndotOmega = n1*omega1 + n2*omega2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotOmega)

            # get tabulated W values for different basis functions np,nq
            Wp = tabWMat[np,kuval,kvval]
            Wq = tabWMat[nq,kuval,kvval]

            if ndim==2
                # Local increment in the location (u,v)
                res += pref*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq

            else
                # add in extra Lval from the action-space volume element (Hamilton et al. 2018, eq 30)
                res += pref*Lval*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq # Local increment in the location (u,v)
            end

        end

        # complete the integration
        res *= deltav
        tabGXi[kuval] = res

    end
    return tabGXi

end



"""
    RunGfuncIsochrone(inputfile)

"""
function RunGfuncIsochrone(inputfile::String,
                           VERBOSE::Int64=1)

    LoadConfiguration(inputfile)

    # Check directory names
    checkdirs = CheckConfigurationDirectories(wmatdir=wmatdir,gfuncdir=gfuncdir)
    if checkdirs < 0
        return 0
    end

    # legendre integration prep
    tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
    tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

    # number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)
    if VERBOSE > 0
        println("CallAResponse.GFuncIsochrone.RunGfuncIsochrone: Considering $nbResVec resonances.")
    end

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)


    Threads.@threads for i = 1:nbResVec
        n1,n2 = tabResVec[1,i],tabResVec[2,i]
        if VERBOSE > 0
            println("CallAResponse.GFuncIsochrone.RunGfuncIsochrone: Starting on ($n1,$n2).")
        end

        # could compute the (u,v) boundaries here (or at least wmin,wmax)
        # compute the frequency scaling factors for this resonance
        ωmin,ωmax = OrbitalElements.FindWminWmaxIsochrone(n1,n2)
        if VERBOSE > 0
            println("CallAResponse.GFuncIsochrone.RunGfuncIsochrone: ωmin=$ωmin,ωmax=$ωmax")
        end

        # for some threading reason, make sure K_u is defined here
        K_u = length(tabwGLquad)

        # loop through once and design a v array for min, max
        vminarr,vmaxarr = zeros(K_u),zeros(K_u)
        for uval = 1:K_u
           vminarr[uval],vmaxarr[uval] = OrbitalElements.FindVminVmaxIsochrone(n1,n2,tabuGLquad[uval])
        end

        # load a value of tabWmat, plus (a,e) values
        filename = wmat_filename(wmatdir,modelname,lharmonic,n1,n2,rb)
        file = h5open(filename,"r")
        Wtab = read(file,"wmat")
        atab = read(file,"amat")
        etab = read(file,"emat")
        nradial,K_u,K_v = size(Wtab)

        # print the size of the found files if the first processor
        if (i==0) & (VERBOSE > 0)
            println("CallAResponse.GFuncIsochrone.RunGfuncIsochrone: Found nradial=$nradial,K_u=$K_u,K_v=$K_v")
        end

        outputfilename = GFuncFilename(gfuncdir,modelname,dfname,lharmonic,n1,n2,K_u,rb)
        if isfile(outputfilename)
            println("CallAResponse.GFuncIsochrone.RunGfuncIsochrone: file already exists for step $i of $nbResVec, ($n1,$n2).")
            continue
        end

        # is having the file open bad?
        # need to loop through all combos of np and nq to make the full matrix.
        h5open(outputfilename, "w") do file

            # loop through all basis function combinations
            for np = 1:nradial
                for nq = 1:nradial

                    if (np==1) & (nq==2) # get a sense of the timing. not the first one, that has compilation time too
                        @time tabGXi = MakeGuIsochrone(ndFdJ,
                                                       n1,n2,np,nq,
                                                       Wtab,atab,etab,
                                                       tabuGLquad,K_v,nradial,
                                                       ωmin,ωmax,
                                                       vminarr,vmaxarr,
                                                       lharmonic,
                                                       ndim=ndim,Omega0=Omega0)
                    else
                        tabGXi = MakeGuIsochrone(ndFdJ,
                                                 n1,n2,np,nq,
                                                 Wtab,atab,etab,
                                                 tabuGLquad,K_v,nradial,
                                                 ωmin,ωmax,
                                                 vminarr,vmaxarr,
                                                 lharmonic,
                                                 ndim=ndim,Omega0=Omega0)
                    end

                    write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)

                end # end nq loop
            end # end np loop

        end # end open file
    end

end
