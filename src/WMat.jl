import OrbitalElements
import AstroBasis
import PerturbPlasma


"""make_wmat(ψ,dψ/dr,d²ψ/dr²,n1,n2,K_u,K_v,lharmonic,basis[,Omega0,NstepsWMat])

@IMPROVE: give basis a type?
@IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
@IMPROVE: parallelise by launching from both -1 and 1?
@IMPROVE: adaptively check for NaN values?

@WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl

@NOTE: AstroBasis provides the Basis_type data type.
"""
function make_wmat(potential::Function,dpotential::Function,ddpotential::Function,n1::Int64,n2::Int64,Kuvals::Matrix{Float64},K_v::Int64,lharmonic::Int64,basis::AstroBasis.Basis_type,Omega0::Float64=1.,NstepsWMat::Int64=20)
    #=
    add rmax as a parameter?
    =#
    K_u = length(Kuvals)

    # compute the frequency scaling factors for this resonance
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    # define beta_c
    beta_c = OrbitalElements.make_betac(dpotential,ddpotential,2000,Omega0)

    # allocate the results matrix
    tabWMat = zeros(basis.nmax,K_u,K_v)
    tabaMat = zeros(K_u,K_v)
    tabeMat = zeros(K_u,K_v)

    # set the matrix step size
    duWMat = (2.0)/(NstepsWMat)

    # start the loop
    for kuval in 1:K_u

        # get the current u value
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

        for kvval in 1:K_v

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e)
            # need to crank the tolerance here, and also check that ecc < 0.
            # put a guard in place for frequency calculations!!
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2)

            # new, iterative brute force procedure
            a1,e1 = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2,1*10^(-12),1)
            maxestep = 0.005
            sma,ecc = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),0)

            # save (a,e) values for later
            tabaMat[kuval,kvval] = sma
            tabeMat[kuval,kvval] = ecc

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            if (rp>ra)
                println("Invalid (rp,ra)=(",rp,",",ra,"). Reversing...check a=",sma," e=",ecc)
                ra,rp = OrbitalElements.rpra_from_ae(sma,ecc)
            end

            # need angular momentum
            Lval = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary by finding Theta(u=-1.)
            gval = OrbitalElements.Theta(potential,dpotential,ddpotential,u,rp,ra,0.02)

            # Uses the Rozier 2019 notation for the mapping to u
            Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5

            # Current location of the radius, r=r(u): isn't this exactly rp?
            rval = Sigma + Delta*OrbitalElements.henon_f(u)

            # the velocity for integration
            dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

            # collect the basis elements (in place!)
            AstroBasis.tabUl!(basis,lharmonic,rval)

            # start the integration loop now that we are initialised
            # at each step, we are performing an RK4-like calculation
            for istep=1:NstepsWMat

                # compute the first prefactor
                pref1 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref1*basis.tabUl[np]
                end

                # update velocities at end of step 1
                k1_1 = duWMat*dt1du
                k2_1 = duWMat*dt2du

                # Step 2
                u += 0.5*duWMat                                                  # Update the time by half a timestep
                rval = Sigma + Delta*OrbitalElements.henon_f(u)                  # Current location of the radius, r=r(u)
                gval = OrbitalElements.Theta(potential,dpotential,ddpotential,u,rp,ra,0.01)
                dt1du, dt2du = omega1*gval, (omega2 - lharmonic/(rval^(2)))*gval # Current value of dtheta1/du and dtheta2/du, always well-posed

                # recompute the basis functions for the changed radius value
                AstroBasis.tabUl!(basis,lharmonic,rval)
                pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_1) + n2*(theta2+0.5*k2_1)) # Common prefactor for all the increments@ATTENTION Depends on the updated (theta1+0.5*k1_1,theta2+0.5*k2_1) !! ATTENTION, to the factor (1.0/3.0) coming from RK4

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref2*basis.tabUl[np]
                end

                # update velocities at end of step 2
                k1_2 = duWMat*dt1du
                k2_2 = duWMat*dt2du

                # Step 3
                # The time, u, is not updated for this step
                # For this step, no need to re-compute the basis elements, as r has not been updated
                pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_2) + n2*(theta2+0.5*k2_2))
                # Common prefactor for all the increments@ATTENTION Depends on the updated (theta1+0.5*k1_2,theta2+0.5*k2_2) !! ATTENTION to the factor (1.0/3.0) coming from RK4

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
                end

                k1_3 = k1_2 # Does not need to be updated
                k2_3 = k2_2 # Does not need to be updated

                # Step 4
                u += 0.5*duWMat # Updating the time by half a timestep: we are now at the next u value
                rval = Sigma + Delta*OrbitalElements.henon_f(u) # Current location of the radius, r=r(u)
                gval = OrbitalElements.Theta(potential,dpotential,ddpotential,u,rp,ra,0.01)
                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval # Current value of dtheta1/du and dtheta2/du, always well-posed

                # updated basis elements for new rval
                AstroBasis.tabUl!(basis,lharmonic,rval)

                pref4 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+k1_3) + n2*(theta2+k2_3)) # Common prefactor for all the increments@ATTENTION Depends on the updated (theta1+k1_3,theta2+k2_3) !! ATTENTION to the factor (1.0/6.0) coming from RK4

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref4*basis.tabUl[np]
                end

                k1_4 = duWMat*dt1du # Current velocity for theta1
                k2_4 = duWMat*dt2du # Current velocity for theta2

                # Update the positions using RK4-like sum
                theta1 += (k1_1 + 2.0*k1_2 + 2.0*k1_3 + k1_4)/(6.0)
                theta2 += (k2_1 + 2.0*k2_2 + 2.0*k2_3 + k2_4)/(6.0)

                # clean nans?

            end # RK4 integration
        end
    end
    return tabWMat,tabaMat,tabeMat
end




"""make_uv_to_ae(ψ,dψ/dr,d²ψ/dr²,n1,n2,K_u,K_v,lharmonic[,Omega0,rb,NstepsWMat])

"""
function make_uv_to_ae(potential::Function,dpotential::Function,ddpotential::Function,n1::Int64,n2::Int64,Kuvals::Matrix{Float64},K_v::Int64,lharmonic::Int64,basis,Omega0::Float64=1.,NstepsWMat::Int64=20)
    #=
    add rmax as a parameter?
    =#
    K_u = length(Kuvals)

    # compute the frequency scaling factors for this resonance
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dpotential,ddpotential,1000.,Omega0)

    # define beta_c
    beta_c = OrbitalElements.make_betac(dpotential,ddpotential,2000,Omega0)

    # allocate the results matrix
    tabaMat = zeros(K_u,K_v)
    tabeMat = zeros(K_u,K_v)

    # set the matrix step size
    duWMat = (2.0)/(NstepsWMat)

    # start the loop
    for kuval in 1:K_u

        # get the current u value
        uval = Kuvals[kuval]

        # get the corresponding v values
        vbound = OrbitalElements.find_vbound(n1,n2,dpotential,ddpotential,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        for kvval in 1:K_v

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dpotential,ddpotential)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e)
            # need to crank the tolerance here, and also check that ecc < 0.
            # put a guard in place for frequency calculations!!
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2)
            a1,e1 = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2,1*10^(-12),1)
            maxestep = 0.005
            sma,ecc = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,omega1,omega2,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),0)


            tabaMat[kuval,kvval] = sma
            tabeMat[kuval,kvval] = ecc
        end
    end
    return tabaMat,tabeMat
end
