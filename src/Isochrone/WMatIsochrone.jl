"""The same as WMat.jl, but Isochrone ψ specific for testing

What's different:
-Use analytic inversion (isochrone_ae_from_omega1omega2)
-Use analytic relationship between alpha and beta for circular orbits

"""

import OrbitalElements
import AstroBasis
import PerturbPlasma


"""MakeWmatIsochrone(ψ,dψ/dr,d²ψ/dr²,n1,n2,K_u,K_v,lharmonic,basis[,Omega0,K_w])

@IMPROVE: give basis a type?
@IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
@IMPROVE: parallelise by launching from both -1 and 1?
@IMPROVE: adaptively check for NaN values?

@WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl
"""
function MakeWmatIsochrone(ψ::Function,dψ::Function,d2ψ::Function,
                           n1::Int64,n2::Int64,
                           Kuvals::Matrix{Float64},
                           K_v::Int64,
                           lharmonic::Int64,
                           basis::AstroBasis.Basis_type,
                           Omega0::Float64=1.,
                           K_w::Int64=20;
                           bc::Float64=1.0,M::Float64=1.0,G::Float64=1.0,EDGE::Float64=0.01)
    #=
    add rmax as a parameter?
    =#
    K_u = length(Kuvals)

    # compute the frequency scaling factors for this resonance
    # not exact, has root finding
    w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dψ,d2ψ,1000.,Omega0)

    # define beta_c: write in for the isochrone
    beta_c = OrbitalElements.analytic_beta_c#OrbitalElements.make_betac(dψ,d2ψ,2000,Omega0)

    # allocate the results matrix
    tabWMat = zeros(basis.nmax,K_u,K_v)
    tabaMat = zeros(K_u,K_v)
    tabeMat = zeros(K_u,K_v)

    # set the matrix step size
    duWMat = (2.0)/(K_w)

    # start the loop
    for kuval in 1:K_u

        # get the current u value
        uval = Kuvals[kuval]

        # get the corresponding v values
        # this requires a zero-finding.
        vbound = OrbitalElements.find_vbound(n1,n2,dψ,d2ψ,1000.,Omega0)
        vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        for kvval in 1:K_v

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)

            # these are exact
            alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dψ,d2ψ)
            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            # convert from omega1,omega2 to (a,e)
            sma,ecc = OrbitalElements.isochrone_ae_from_omega1omega2(omega1,omega2,bc,M,G)

            # save (a,e) values for later
            tabaMat[kuval,kvval] = sma
            tabeMat[kuval,kvval] = ecc

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need angular momentum
            Lval = OrbitalElements.isochrone_L_from_rpra(rp,ra,bc,M,G)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary by finding ThetaRpRa(u=-1.)
            gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)

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
            for istep=1:K_w

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
                # Update the time by half a timestep
                u += 0.5*duWMat

                # Current location of the radius, r=r(u)
                rval = Sigma + Delta*OrbitalElements.henon_f(u)

                # current value of Theta
                gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)

                # Current value of dtheta1/du and dtheta2/du
                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

                # recompute the basis functions for the changed radius value
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # Common prefactor for all the increments@ATTENTION Depends on the updated (theta1+0.5*k1_1,theta2+0.5*k2_1)
                # the factor (1.0/3.0) coming from RK4
                pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_1) + n2*(theta2+0.5*k2_1))

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
                gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
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


"""
    RunWmatIsochrone(inputfile)

"""
function RunWmatIsochrone(inputfile::String)

    # load model parameters
    include(inputfile)

    # check directory before proceeding (save time if not.)
    if !(isdir(wmatdir))
        error("WMat.jl:: wmatdir not found")
    end

    # bases prep.
    AstroBasis.fill_prefactors!(basis)
    bases=[deepcopy(basis) for k=1:Threads.nthreads()]

    # Legendre integration prep.
    tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
    tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

    # number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # print the length of the list of resonance vectors
    println("WMat.jl: Number of resonances to compute: $nbResVec")

    Threads.@threads for i = 1:nbResVec
        k = Threads.threadid()
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        println("WMat.jl: Computing W for the ($n1,$n2) resonance.")

        # currently defaulting to timed version:
        # could make this a flag (timing optional)
        @time tabWMat,tabaMat,tabeMat = MakeWmatIsochrone(ψ,dψ,d2ψ,n1,n2,tabuGLquad,K_v,lharmonic,bases[k],Omega0,K_w,bc=bc,M=M,G=G)

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(wmat_filename(wmatdir,modelname,lharmonic,n1,n2,rb), "w") do file
            write(file, "wmat",tabWMat)
            write(file, "amat",tabaMat)
            write(file, "emat",tabeMat)
        end

    end

end
