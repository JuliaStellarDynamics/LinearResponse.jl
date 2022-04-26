
# we need to align the G to be a simple function of u

"""
wrap the ndFdJ function you want
"""
function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    return ROIndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)
    #return ISOndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG)
end


"""
# FROM MEAN.JL
# Function that returns the value of n.dF/dJ
# where the DF follows the normalisation convention int dxdv F = Mtot
# Arguments are:
# + (n1,n2)
# + (E,L)
# + ndotOmega = n1*Omega1 + n2*Omega2
# ATTENTION, specific to an isotropic DF
"""
function ISOndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG) # Current value of dF/E. ATTENTION, the DF is assumed to be isotropic
    res = ndotOmega*dFdE # Current value of n.dF/dJ. ATTENTION, the DF is assumed to be isotropic
    #####
    return res
end

"""
# FROM ROI.JL
# Function that computes n.dF/dJ for the ROI DF
# ATTENTION, we put n.dF/dJ != 0 only for 0 < Q < 1
# @IMPROVE -- ideally we should have rather reduced
# the (alpha,beta) domain of integration
"""
function ROIndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)

    Q = OrbitalElements.isochrone_Q_ROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochrone_Saha_dDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochrone_dQdE_ROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochrone_dQdL_ROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res
end


function make_Gu(potential::Function,dpotential::Function,ddpotential::Function,n1::Int64,n2::Int64,np::Int64,nq::Int64,tabWMat::Array{Float64},tabaMat::Array{Float64},tabeMat::Array{Float64},Kuvals::Matrix{Float64},K_v::Int64,nradial::Int64,lharmonic::Int64,Omega0::Float64=1.,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)


    # get prefactors
    CMatrix = getCMatrix(lharmonic)

    # get basic parameters
    K_u     = length(Kuvals)

    # set up a blank array
    tabGXi  = zeros(K_u)

    # compute the frequency scaling factors for this resonance
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    # define beta_c
    beta_c = OrbitalElements.make_betac(dpotential,ddpotential,2000,Omega0)

    # Overall prefactor of the integrand. ATTENTION, to the minus sign.
    pref = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*2+1.0)

    for kuval in 1:K_u

        uval = Kuvals[kuval]

        # get the corresponding v values
        #omega1_func(x) = Omega1_circular(dpotential,ddpotential,x)
        # m is the extremum location of omegafunc: get from find_wmin_wmax?
        # OmegaC is the central Omega value, which is Omega0 (I think?)
        #vbound = omega1_func(m)/Omega0
        #extreme(x) = n1*OrbitalElements.Omega1_circular(dpotential,ddpotential,x) + #n2*OrbitalElements.Omega2_circular(dpotential,x)
        #m = OrbitalElements.extremise_function(extreme,24,0.,1000.,false)
        #vbound = OrbitalElements.Omega1_circular(dpotential,ddpotential,m)/Omega0
        #vbound =
        vbound = OrbitalElements.find_vbound(n1,n2,dpotential,ddpotential,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        res = 0.0 # Initialising the result

        for kvval in 1:K_v
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e)
            sma,ecc  = tabaMat[kuval,kvval],tabeMat[kuval,kvval]

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need (E,L)
            Lval = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
            Eval = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

            # compute Jacobians
            Jacalphabeta = OrbitalElements.Jacalphabeta_to_uv(n1,n2,w_min,w_max,0.5) #(alpha,beta) -> (u,v)
            JacEL        = OrbitalElements.isochrone_JacEL_to_alphabeta(alpha,beta)          #(E,L) -> (alpha,beta) # ISOCHRONE SPECIFIC!!
            JacJ         = (1/omega1)                                #(J) -> (E,L)
            dimensionl   = (1/Omega0)                                # remove dimensionality


            # get the resonance vector
            ndotOmega = n1*omega1 + n2*omega2

            # compute dF/dJ: this can/should be redefined
            valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotOmega)

            # get tabulated W values for different basis functions np,nq
            Wp = tabWMat[np,kuval,kvval]
            Wq = tabWMat[nq,kuval,kvval]

            # do a nan check?
            nancheck = true
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
