

"""MakeWmatUV(ψ,dψ,d2ψ,d3ψ,n1,n2,K_u,K_v,lharmonic,basis[,Ω0,K_w])

@IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
@IMPROVE: parallelise by launching from both -1 and 1?
@IMPROVE: adaptively check for NaN values?

@WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl

@NOTE: AstroBasis provides the Basis_type data type.
"""
function MakeWmatUV(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    n1::Int64,n2::Int64,
                    Kuvals::Matrix{Float64},
                    K_v::Int64,
                    lharmonic::Int64,
                    basis::AstroBasis.Basis_type;
                    Ω0::Float64=1.,
                    rmin::Float64=1.0e-6,
                    rmax::Float64=1.0e4,
                    K_w::Int64=50,
                    EDGE::Float64=0.01,
                    VERBOSE::Int64=0,
                    NINT::Int64=32)

    # get the number of u samples from the input vector of u vals
    K_u = length(Kuvals)

    # Frequency cuts associated to [rmin,rmax]
    # @IMPROVE: compute them once (independant of n1,n2) and function argument ?
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω0=Ω0)
    # compute the frequency scaling factors for this resonance
    ωmin,ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,Ω0=Ω0,rmin=rmin,rmax=rmax)

    # define a function for beta_c
    βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω0,rmin=rmin,rmax=rmax)

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

        if VERBOSE > 1
            println("\nCallAResponse.WMat.MakeWMat: on step $kuval of $K_u: u=$uval.")
        end

        # get the corresponding v boundary values
        vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω0=Ω0,rmin=rmin,rmax=rmax)

        # determine the step size in v
        deltav = (vmax - vmin)/(K_v)

        for kvval in 1:K_v

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            # big step: convert input (u,v) to (rp,ra)
            # now we need (rp,ra) that corresponds to (u,v)
            α,β = OrbitalElements.AlphaBetaFromUV(uval,vval,n1,n2,ωmin,ωmax)

            Ω1,Ω2 = α*Ω0,α*β*Ω0

            if VERBOSE > 2
                print("v=$kvval,o1=$Ω1,o2=$Ω2;")
            end

            # convert from Ω1,Ω2 to (a,e)
            # need to crank the tolerance here, and also check that ecc < 0.
            # put a guard in place for frequency calculations!!
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,Ω1,Ω2)

            # new, iterative brute force procedure
            #a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,Ω1,Ω2,1*10^(-12),1)
            #maxestep = 0.005
            #sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,Ω1,Ω2,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),0)
            a,e = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,NINT=NINT,EDGE=EDGE,VERBOSE=VERBOSE)

            #if VERBOSE > 2
            #    print("a=$sma,ecc=$ecc")
            #end

            # save (a,e) values for later
            tabaMat[kuval,kvval] = a
            tabeMat[kuval,kvval] = e

            # compute the Jacobian (E,L)->(alpha,beta) here. a little more expensive, but savings in the long run
            #tabJMat[kuval,kvval] = OrbitalElements.JacELToAlphaBetaAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,Ω0=Ω0)
            tabJMat[kuval,kvval] = OrbitalElements.JacELToAlphaBetaAE(a,e,ψ,dψ,d2ψ,Ω0)

            # need angular momentum
            Lval = OrbitalElements.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary by finding ThetaRpRa(u=-1.)
            #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
            gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

            # Current location of the radius, r=r(u): isn't this exactly rp?
            rval = OrbitalElements.ru(u,a,e)

            # the velocity for integration
            dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

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
                rval = OrbitalElements.ru(u,a,e)

                # get the corresponding value of Theta(u): could use (a,e) function here as well
                #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
                gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

                # Current value of dtheta1/du and dtheta2/du
                dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

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
                # Current location of the radius, r=r(u)
                rval = OrbitalElements.ru(u,a,e)

                # current value of dtheta1/du and dtheta2/du
                #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
                gval = OrbitalElements.ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

                dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

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
    RunWmat()

"""
function RunWmat(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                 wmatdir::String,
                 K_u::Int64,K_v::Int64,K_w::Int64,
                 basis::AstroBasis.Basis_type,
                 lharmonic::Int64,
                 n1max::Int64,
                 nradial::Int64,
                 Ω0::Float64,
                 modelname::String,
                 rb::Float64,
                 rmin::Float64,rmax::Float64;
                 VERBOSE::Int64=0)

    # load model parameters
    #include(inputfile)

    # check wmat directory before proceeding (save time if not.)
    checkdirs = CheckConfigurationDirectories(wmatdir=wmatdir)
    if checkdirs < 0
        return 0
    end

    # get basis parameters
    ndim = basis.dimension
    nradialmax = basis.nmax

    # check if we can cover the specified radial orders
    if nradialmax > nradial
        println("CallAResponse.WMat.RunWmat: the input basis does not have sufficient nradial ($nradialmax) for the requested value ($nradial).")
        return 0
    end

    # bases prep.
    AstroBasis.fill_prefactors!(basis)
    bases=[deepcopy(basis) for k=1:Threads.nthreads()]

    # Legendre integration prep.
    tabuGLquadtmp,tabwGLquad = FiniteHilbertTransform.tabuwGLquad(K_u)
    tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

    # number of resonance vectors
    nbResVec = get_nbResVec(lharmonic,n1max,ndim)

    # fill in the array of resonance vectors (n1,n2)
    tabResVec = maketabResVec(lharmonic,n1max,ndim)

    # print the length of the list of resonance vectors
    if VERBOSE>0
        println("CallAResponse.WMat.RunWmat: Number of resonances to compute: $nbResVec")
    end

    Threads.@threads for i = 1:nbResVec
        k = Threads.threadid()
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        if VERBOSE>0
            println("CallAResponse.WMat.RunWmat: Computing W for the ($n1,$n2) resonance.")
        end

        if isfile(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,nradial,rb,K_u,K_v,K_w))
            if VERBOSE > 0
                println("CallAResponse.WMatIsochrone.RunWmatIsochrone: ($n1,$n2) resonanance WMat file already exists.")
            end
            continue
        end

        # compute the W matrices in UV space: timing optional
        if VERBOSE>1
            @time tabWMat,tabaMat,tabeMat,tabJMat = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,
                                                               n1,n2,
                                                               tabuGLquad,K_v,lharmonic,bases[k],
                                                               Ω0=Ω0,rmin=rmin,rmax=rmax,K_w=K_w,VERBOSE=VERBOSE)
        else
            tabWMat,tabaMat,tabeMat,tabJMat = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,
                                                         n1,n2,
                                                         tabuGLquad,K_v,lharmonic,bases[k],
                                                         Ω0=Ω0,rmin=rmin,rmax=rmax,K_w=K_w,VERBOSE=VERBOSE)
        end


        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,nradial,rb,K_u,K_v,K_w), "w") do file
            write(file, "wmat",tabWMat)
            write(file, "amat",tabaMat)
            write(file, "emat",tabeMat)
            write(file, "jELABmat",tabJMat)
        end

    end

end
