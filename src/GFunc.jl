


"""
function to compute G(u)

@ATTENTION, the dimensionality (e.g. 2d vs 3d) is now encoded in 'ndim'.

"""
function makeGu(potential::Function,dpotential::Function,ddpotential::Function,
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
        # 2d prefactor, see Fouvry et al. 2015
        pref = (2*pi)^2
    else
        # 3d prefactor, see Hamilton et al. 2018
        CMatrix = getCMatrix(lharmonic)
        pref    = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*lharmonic+1.0)
    end

    # get basic parameters
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    # compute the frequency scaling factors for this resonance
    ωmin,ωmax = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    # define beta_c
    beta_c = OrbitalElements.make_betac(dpotential,ddpotential,2000,Omega0)

    for kuval in 1:K_u

        uval = Kuvals[kuval]

        vbound = OrbitalElements.find_vbound(n1,n2,dpotential,ddpotential,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,ωmin,ωmax,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        res = 0.0 # Initialising the result

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            #alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential,1000.,Omega0)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,ωmin,ωmax)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e): using a tabled value
            sma,ecc  = tabaMat[kuval,kvval],tabeMat[kuval,kvval]

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need (E,L)
            Lval = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
            Eval = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

            # compute Jacobians
            Jacalphabeta = OrbitalElements.Jacalphabeta_to_uv(n1,n2,ωmin,ωmax,vval) #(alpha,beta) -> (u,v). owing to the remapping of omega, this has an extra 2/(ωmax-ωmin)
            JacEL        = OrbitalElements.JacEL_to_alphabeta(alpha,beta)          #(E,L) -> (alpha,beta)
            JacJ         = (1/omega1)                                #(J) -> (E,L)
            dimensionl   = (1/Omega0)                                # remove dimensionality from omega mapping


            # get the resonance vector
            ndotOmega = n1*omega1 + n2*omega2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotOmega)

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
                    println(Jacalphabeta," ",JacEL," ",pref," ",(Lval/omega1)," ",valndFdJ," ",Wp," ",Wq)
                end
            end
            =#

            if ndim==2
                res += pref*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq # Local increment in the location (u,v)

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
