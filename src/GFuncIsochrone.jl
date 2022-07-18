


"""
function to compute G(u)

@ATTENTION, the dimensionality (e.g. 2d vs 3d) is now encoded in 'pref'.

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
                 lharmonic::Int64,
                 pref::Float64;
                 Omega0::Float64=1.,bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    # get basic parameters
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    # compute the frequency scaling factors for this resonance
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    # define beta_c
    beta_c = OrbitalElements.make_betac(dpotential,ddpotential,2000,Omega0)

    for kuval in 1:K_u

        uval = Kuvals[kuval]

        vbound = OrbitalElements.find_vbound(n1,n2,dpotential,ddpotential,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        res = 0.0 # Initialising the result

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential,1000.,Omega0)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e)
            #sma,ecc  = tabaMat[kuval,kvval],tabeMat[kuval,kvval]
            sma,ecc = OrbitalElements.isochrone_ae_from_omega1omega2(omega1,omega2,bc,M,G)

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need (E,L)
            Lval = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
            Eval = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

            # compute Jacobians
            Jacalphabeta = OrbitalElements.Jacalphabeta_to_uv(n1,n2,w_min,w_max,vval) #(alpha,beta) -> (u,v)
            JacEL        = OrbitalElements.JacEL_to_alphabeta(alpha,beta)          #(E,L) -> (alpha,beta)
            JacJ         = (1/omega1)                                #(J) -> (E,L)
            dimensionl   = (1/Omega0)                                # remove dimensionality


            # get the resonance vector
            ndotOmega = n1*omega1 + n2*omega2

            # compute dF/dJ: call out for value
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotOmega)

            # get tabulated W values for different basis functions np,nq
            Wp = tabWMat[np,kuval,kvval]
            Wq = tabWMat[nq,kuval,kvval]

            # do a nan check?
            nancheck = false
            if (nancheck)
                tmp = pref*lharmonic*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq

                if isnan(tmp)
                    println(Jacalphabeta," ",JacEL," ",pref," ",(Lval/omega1)," ",valndFdJ," ",Wp," ",Wq)
                end
            end

            res += pref*Lval*(dimensionl*Jacalphabeta*JacEL*JacJ*valndFdJ)*Wp*Wq # Local increment in the location (u,v)

        end

        # complete the integration
        res *= deltav
        tabGXi[kuval] = res

    end
    return tabGXi

end
