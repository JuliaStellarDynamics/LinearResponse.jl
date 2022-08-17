import OrbitalElements
import AstroBasis
import PerturbPlasma
using HDF5


"""MakeWmat(ψ,dψ,d2ψ,d3ψ,n1,n2,K_u,K_v,lharmonic,basis[,Omega0,K_w])

@IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
@IMPROVE: parallelise by launching from both -1 and 1?
@IMPROVE: adaptively check for NaN values?

@WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl

@NOTE: AstroBasis provides the Basis_type data type.
"""
function MakeWmat(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                  n1::Int64,n2::Int64,
                  Kuvals::Matrix{Float64},
                  K_v::Int64,
                  lharmonic::Int64,
                  basis::AstroBasis.Basis_type;
                  Omega0::Float64=1.,
                  K_w::Int64=50,
                  EDGE::Float64=0.01,
                  verbose::Int64=0)
    #=
    add rmax as a parameter?
    =#

    # get the number of u samples from the input vector of u vals
    K_u = length(Kuvals)

    # compute the frequency scaling factors for this resonance
    # @IMPROVE: need a maximum radius for the root finding; here set to be 1000., but if the cluster was extremely extended, this could break
    ωmin,ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,10000.,Omega0)

    # define beta_c, empirically.
    # @IMPROVE: 2000 is the number of sample points; this also is hard-coded to sample between log R [-5,5], which would be adaptive
    # beta_c = OrbitalElements.make_betac(dψ,d2ψ,2000,Omega0)
    beta_c(alpha_c::Float64)::Float64 = OrbitalElements.beta_circ(alpha_c,dψ,d2ψ,Omega0,rmax=1.0e6)

    # allocate the results matrices
    tabWMat = zeros(basis.nmax,K_u,K_v)
    tabaMat = zeros(K_u,K_v)
    tabeMat = zeros(K_u,K_v)
    tabJMat = zeros(K_u,K_v)

    # set the matrix step size
    duWMat = (2.0)/(K_w)

    # start the loop
    for kuval in 1:K_u

        # get the current u value
        uval = Kuvals[kuval]

        if verbose > 1
            println("\nCallAResponse.WMat.MakeWMat: on step $kuval of $K_u: u=$uval.")
        end

        # get the corresponding v boundary values
        vbound = OrbitalElements.FindVbound(n1,n2,dψ,d2ψ,10000.,Omega0)
        vmin,vmax = OrbitalElements.FindVminVmax(uval,ωmin,ωmax,n1,n2,vbound,beta_c)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        for kvval in 1:K_v

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)



            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            alpha,beta = OrbitalElements.AlphaBetaFromUV(uval,vval,n1,n2,ωmin,ωmax)

            omega1,omega2 = alpha*Omega0,alpha*beta*Omega0

            if verbose > 2
                print("v=$kvval,o1=$omega1,o2=$omega2;")
            end

            # convert from omega1,omega2 to (a,e)
            # need to crank the tolerance here, and also check that ecc < 0.
            # put a guard in place for frequency calculations!!
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,omega1,omega2)

            # new, iterative brute force procedure
            #a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,omega1,omega2,1*10^(-12),1)
            #maxestep = 0.005
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,omega1,omega2,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),0)
            sma,ecc = OrbitalElements.AEFromOmega1Omega2Brute(omega1,omega2,ψ,dψ,d2ψ,d3ψ,NINT=32,EDGE=EDGE,verbose=0)

            #if verbose > 2
            #    print("a=$sma,ecc=$ecc")
            #end

            # save (a,e) values for later
            tabaMat[kuval,kvval] = sma
            tabeMat[kuval,kvval] = ecc

            # compute the Jacobian (E,L)->(alpha,beta) here. a little more expensive, but savings in the long run
            #tabJMat[kuval,kvval] = OrbitalElements.JacELToAlphaBetaAE(ψ,dψ,d2ψ,d3ψ,d4ψ,sma,ecc,Omega0=Omega0)
            tabJMat[kuval,kvval] = OrbitalElements.JacELToAlphaBetaAE(sma,ecc,ψ,dψ,d2ψ,Omega0)

            # get (rp,ra)
            rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)

            # need angular momentum
            Lval = OrbitalElements.LFromRpRa(ψ,dψ,d2ψ,rp,ra)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary by finding ThetaRpRa(u=-1.)
            #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
            gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,sma,ecc,EDGE=EDGE)

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

                # Begin RK4 Step 2
                # Update the time by half a timestep
                u += 0.5*duWMat

                # Current location of the radius, r=r(u)
                rval = Sigma + Delta*OrbitalElements.henon_f(u)

                # get the corresponding value of Theta(u): could use (a,e) function here as well
                #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
                gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,sma,ecc,EDGE=EDGE)

                # Current value of dtheta1/du and dtheta2/du
                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

                # recompute the basis functions for the changed radius value
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # common prefactor for all the increments
                # depends on the updated (theta1+0.5*k1_1,theta2+0.5*k2_1)
                # the factor (1.0/3.0) comes from RK4
                pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_1) + n2*(theta2+0.5*k2_1))

                # Loop over the radial indices to sum basis contributions
                #for np=1:basis.nmax
                #    tabWMat[np,kuval,kvval] += pref2*basis.tabUl[np]
                #end

                # update velocities at end of step 2
                k1_2 = duWMat*dt1du
                k2_2 = duWMat*dt2du

                # Begin step 3 of RK4
                # The time, u, is not updated for this step
                # For this step, no need to re-compute the basis elements, as r has not been updated
                # Common prefactor for all the increments
                # Depends on the updated (theta1+0.5*k1_2,theta2+0.5*k2_2)
                # the factor (1.0/3.0) comes from RK4
                pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_2) + n2*(theta2+0.5*k2_2))

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    #tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
                    tabWMat[np,kuval,kvval] += (pref2+pref3)*basis.tabUl[np]
                end

                k1_3 = k1_2 # Does not need to be updated
                k2_3 = k2_2 # Does not need to be updated

                # Begin step 4 of RK4
                u += 0.5*duWMat # Updating the time by half a timestep: we are now at the next u value
                rval = Sigma + Delta*OrbitalElements.henon_f(u) # Current location of the radius, r=r(u)

                # current value of dtheta1/du and dtheta2/du
                #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
                gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,sma,ecc,EDGE=EDGE)

                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

                # updated basis elements for new rval
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # Common prefactor for all the increments
                # Depends on the updated (theta1+k1_3,theta2+k2_3)
                # The factor (1.0/6.0) comes from RK4
                pref4 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+k1_3) + n2*(theta2+k2_3))

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref4*basis.tabUl[np]
                end

                k1_4 = duWMat*dt1du # Current velocity for theta1
                k2_4 = duWMat*dt2du # Current velocity for theta2

                # Update the positions using RK4-like sum
                theta1 += (k1_1 + 2.0*k1_2 + 2.0*k1_3 + k1_4)/(6.0)
                theta2 += (k2_1 + 2.0*k2_2 + 2.0*k2_3 + k2_4)/(6.0)

                # clean or check nans?

            end # RK4 integration
        end
    end
    return tabWMat,tabaMat,tabeMat,tabJMat
end


"""
    RunWmat(inputfile)

"""
function RunWmat(inputfile::String)

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
        @time tabWMat,tabaMat,tabeMat,tabJMat = MakeWmat(ψ,dψ,d2ψ,d3ψ,n1,n2,tabuGLquad,K_v,lharmonic,bases[k],Omega0=Omega0,K_w=K_w,verbose=2)

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(wmat_filename(wmatdir,modelname,lharmonic,n1,n2,rb), "w") do file
            write(file, "wmat",tabWMat)
            write(file, "amat",tabaMat)
            write(file, "emat",tabeMat)
            write(file, "jELABmat",tabJMat)
        end

    end

end
