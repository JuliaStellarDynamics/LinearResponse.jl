


"""MakeGuIsochrone
function to compute G(u), isochrone specific.

"""
function MakeGuIsochrone(potential::Function,dpotential::Function,ddpotential::Function,
                         ndFdJ::Function,
                         n1::Int64,n2::Int64,
                         np::Int64,nq::Int64,
                         tabWMat::Array{Float64},
                         tabaMat::Array{Float64},
                         tabeMat::Array{Float64},
                         Kuvals::Matrix{Float64},
                         K_v::Int64,nradial::Int64,
                         lharmonic::Int64;
                         ndim::Int64,
                         Omega0::Float64=1.,bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    # calculate the prefactor based on the dimensionality (defaults to 3d)
    if ndim==2
        pref = (2*pi)^2
    else
        CMatrix = getCMatrix(lharmonic)
        pref    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*lharmonic+1.0)
    end

    # get basic parameters
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    # compute the frequency scaling factors for this resonance
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    for kuval in 1:K_u

        uval = Kuvals[kuval]

        # get the v integration boundaries. these are not perfectly exact (requires a zero-finding), but all other items going into the calculation are.
        vbound = OrbitalElements.find_vbound(n1,n2,dpotential,ddpotential,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,OrbitalElements.analytic_beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        res = 0.0 # Initialising the result

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential,1000.,Omega0)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e) using isochrone exact version
            sma,ecc = OrbitalElements.isochrone_ae_from_omega1omega2(omega1,omega2,bc,M,G)

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need (E,L), use isochrone exact version
            Eval,Lval = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)

            # compute Jacobians
            #(alpha,beta) -> (u,v)
            Jacalphabeta = OrbitalElements.JacalphabetaToUV(n1,n2,w_min,w_max,vval)

            #(E,L) -> (alpha,beta): Isochrone analytic
            JacEL        = OrbitalElements.isochrone_JacEL_to_alphabeta(alpha,beta,bc,M,G)

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
function RunGfuncIsochrone(inputfile::String)

    include(inputfile)

    #####
    # Check directories names
    #####
    if !(isdir(wmatdir) && isdir(gfuncdir))
        error("GFunc.jl: wmatdir or gfuncdir not found ")
    end

    #####
    # Legendre integration prep.
    #####
    tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
    tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

    #####
    # Construct the table of needed resonance vectors
    #####

    # Number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)
    tabResVec = maketabResVec(lharmonic,n1max,ndim) # Filling in the array of resonance vectors (n1,n2)

    println("GFunc.jl: Considering $nbResVec resonances.")

    Threads.@threads for i = 1:nbResVec
        n1,n2 = tabResVec[1,i],tabResVec[2,i]
        println("Gfunc.jl: Starting on ($n1,$n2).")

        # load a value of tabWmat, plus (a,e) values
        filename = wmat_filename(wmatdir,modelname,lharmonic,n1,n2,rb)
        file = h5open(filename,"r")
        Wtab = read(file,"wmat")
        atab = read(file,"amat")
        etab = read(file,"emat")
        nradial,K_u,K_v = size(Wtab)

        # print the size of the found files if the first processor
        if i==0
            println("GFunc.jl: Found nradial=$nradial,K_u=$K_u,K_v=$K_v")
        end

        # need to loop through all combos of np and nq to make the full matrix.
        h5open(gfunc_filename(gfuncdir,modelname,dfname,lharmonic,n1,n2,K_u), "w") do file

            # loop through all basis function combinations
            for np = 1:nradial
                for nq = 1:nradial

                    tabGXi = MakeGuIsochrone(potential,dpotential,ddpotential,ndFdJ,n1,n2,np,nq,Wtab,atab,etab,tabuGLquad,K_v,nradial,lharmonic,ndim=ndim,Omega0=Omega0)
                    sumG = sum(tabGXi)
                    if (np>-100) & (nq>-100)
                        if isnan(sumG)
                            println("NaN for n1=$n1, n2=$n2.")
                        else
                            #println("np=$np, nq=$nq, sumG=$sumG.")
                        end
                    end
                    write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)
                end
            end
        end
    end

end
