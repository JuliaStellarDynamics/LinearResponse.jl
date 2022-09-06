

"""
    MakeWmatUV(ψ,dψ,d2ψ,d3ψ,n1,n2,tabu,Kv,lharmonic,basis[,Ω₀,rmin,rmax,K_w])

@TO DESCRIBE
"""
function MakeWmatUV(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    n1::Int64,n2::Int64,
                    tabu::Array{Float64},
                    Kv::Int64,
                    lharmonic::Int64,
                    basis::AstroBasis.Basis_type;
                    Ω₀::Float64=defaultΩ₀,
                    rmin::Float64=defaultrmin,
                    rmax::Float64=defaultrmax,
                    Kw::Int64=50,
                    EDGE::Float64=0.01,
                    VERBOSE::Int64=0,
                    NINT::Int64=32)

    """
    @IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
    @IMPROVE: parallelise by launching from both -1 and 1?
    @IMPROVE: adaptively check for NaN values?

    @WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl

    @NOTE: AstroBasis provides the Basis_type data type.
    """

    # get the number of u samples from the input vector of u vals
    Ku = length(tabu)

    # Frequency cuts associated to [rmin,rmax]
    # @IMPROVE: compute them once (independant of n1,n2) and function argument ?
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)
    # compute the frequency scaling factors for this resonance
    ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

    ωminmax = [ωmin,ωmax]

    # define a function for beta_c
    βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin=rmin,rmax=rmax)

    # allocate the results matrices
    tabWMat = zeros(basis.nmax,Ku,Kv)
    tabΩ1Ω2Mat = zeros(Ku,Kv,2)
    tabAEMat = zeros(Ku,Kv,2)
    tabJMat = zeros(Ku,Kv)
    tabvminmax = zeros(Ku,2)

    # set the matrix step size
    duWMat = (2.0)/(Kw)

    # start the loop
    for kuval in 1:Ku

        # get the current u value
        uval = tabu[kuval]

        if VERBOSE > 2
            println("\nCallAResponse.WMat.MakeWMat: on step $kuval of $Ku: u=$uval.")
        end

        # get the corresponding v boundary values
        vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

        # saving them
        tabvminmax[kuval,1], tabvminmax[kuval,2] = vmin, vmax

        # determine the step size in v
        deltav = (vmax - vmin)/(Kv)

        for kvval in 1:Kv

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            ####
            # (u,v) -> (a,e)
            ####
            # (u,v) -> (α,β)
            α,β = OrbitalElements.αβFromUV(uval,vval,n1,n2,ωmin,ωmax)
            # (α,β) -> (Ω1,Ω2)
            Ω1,Ω2 = α*Ω₀,α*β*Ω₀
            # (Ω1,Ω2) -> (a,e)
            a,e = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,NINT=NINT,EDGE=EDGE,VERBOSE=VERBOSE)

            if VERBOSE > 2
                print("v=$kvval,o1=$Ω1,o2=$Ω2;")
            end

            # save (Ω1,Ω2) values for later
            tabΩ1Ω2Mat[kuval,kvval,1], tabΩ1Ω2Mat[kuval,kvval,2] = Ω1, Ω2
            # save (a,e) values for later
            tabAEMat[kuval,kvval,1], tabAEMat[kuval,kvval,2] = a, e

            # compute the Jacobian (E,L)->(alpha,beta) here. a little more expensive, but savings in the long run
            tabJMat[kuval,kvval] = OrbitalElements.JacELToαβAE(a,e,ψ,dψ,d2ψ,Ω₀)

            # need angular momentum
            Lval = OrbitalElements.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary
            gval = OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

            # Current location of the radius, r=r(u)
            rval = OrbitalElements.ru(u,a,e)

            # the velocity for integration
            dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

            # collect the basis elements (in place!)
            AstroBasis.tabUl!(basis,lharmonic,rval)

            # start the integration loop now that we are initialised
            # at each step, we are performing an RK4-like calculation
            for istep=1:Kw

                ####
                # RK4 Step 1
                ####
                # compute the first prefactor
                pref1 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref1*basis.tabUl[np]
                end

                # update velocities at end of step 1
                k1_1 = duWMat*dt1du
                k2_1 = duWMat*dt2du

                ####
                # RK4 Step 2
                ####
                # Update the time by half a timestep
                u += 0.5*duWMat

                # Current location of the radius, r=r(u)
                rval = OrbitalElements.ru(u,a,e)

                # get the corresponding value of Theta(u): could use (a,e) function here as well
                gval = OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

                # Current value of dtheta1/du and dtheta2/du
                dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

                # recompute the basis functions for the changed radius value
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # common prefactor for all the increments
                # depends on the updated (theta1+0.5*k1_1,theta2+0.5*k2_1)
                # the factor (1.0/3.0) comes from RK4
                pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_1) + n2*(theta2+0.5*k2_1))

                # update velocities at end of step 2
                k1_2 = duWMat*dt1du
                k2_2 = duWMat*dt2du

                ####
                # RK4 step 3
                ####
                # The time, u, is not updated for this step
                # For this step, no need to re-compute the basis elements, as r has not been updated
                # Common prefactor for all the increments
                # Depends on the updated (theta1+0.5*k1_2,theta2+0.5*k2_2)
                # the factor (1.0/3.0) comes from RK4
                pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k1_2) + n2*(theta2+0.5*k2_2))

                # Loop over the radial indices to sum basis contributions
                # Contribution of steps 2 and 3 together
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += (pref2+pref3)*basis.tabUl[np]
                end

                k1_3 = k1_2 # Does not need to be updated
                k2_3 = k2_2 # Does not need to be updated

                ####
                # RK4 step 4
                ####
                u += 0.5*duWMat # Updating the time by half a timestep: we are now at the next u value
                # Current location of the radius, r=r(u)
                rval = OrbitalElements.ru(u,a,e)

                # current value of dtheta1/du and dtheta2/du
                #gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,u,rp,ra,EDGE=EDGE)
                gval = OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

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
    return tabWMat,tabΩ1Ω2Mat,tabAEMat,tabJMat,tabvminmax,ωminmax
end


"""
    RunWmat(ψ,dψ,d2ψ,d3ψ,wmatdir,Ku,Kv,Kw,basis,lharmonic,n1max,nradial,Ω₀,modelname,rb,rmin,rmax[,VERBOSE])

@TO DESCRIBE
"""
function RunWmat(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                 wmatdir::String,
                 FHT::FiniteHilbertTransform.FHTtype,
                 Kv::Int64,Kw::Int64,
                 basis::AstroBasis.Basis_type,
                 lharmonic::Int64,
                 n1max::Int64,
                 Ω₀::Float64,
                 modelname::String,
                 rmin::Float64,rmax::Float64;
                 VERBOSE::Int64=0,
                 OVERWRITE::Bool=false)

    # check wmat directory before proceeding (save time if not.)
    checkdirs = CheckConfigurationDirectories(wmatdir=wmatdir)
    if checkdirs < 0
        return 0
    end

    # get basis parameters
    ndim, nradial, rb = basis.dimension, basis.nmax, basis.rb

    # bases prep.
    AstroBasis.fill_prefactors!(basis)
    bases=[deepcopy(basis) for k=1:Threads.nthreads()]

    # Integration points
    tabu, Ku = FHT.tabu, FHT.Ku

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

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

        # If it has been already computed
        if isfile(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw))
            file = h5open(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw), "r")
            oldnradial = read(file,"nradial")
            if (OVERWRITE == false) && (nradial <= oldnradial)
                if VERBOSE > 0
                    println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with higher nradial: no computation.")
                end
                continue
            else
                if VERBOSE > 0
                    println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with lower nradial: recomputing and overwritting.")
                end
            end
            close(file)
        end

        # compute the W matrices in UV space: timing optional
        if VERBOSE>1
            @time tabWMat,tabΩ1Ω2Mat,tabAEMat,tabJMat,tabvminmax,ωminmax = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,
                                                                                    n1,n2,
                                                                                    tabu,Kv,lharmonic,bases[k],
                                                                                    Ω₀=Ω₀,rmin=rmin,rmax=rmax,Kw=Kw,VERBOSE=VERBOSE)
        else
            tabWMat,tabΩ1Ω2Mat,tabAEMat,tabJMat,tabvminmax,ωminmax = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,
                                                                            n1,n2,
                                                                            tabu,Kv,lharmonic,bases[k],
                                                                            Ω₀=Ω₀,rmin=rmin,rmax=rmax,Kw=Kw,VERBOSE=VERBOSE)
        end


        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw), "w") do file
            write(file, "nradial",nradial)
            write(file, "wmat",tabWMat)
            write(file, "Omgmat",tabΩ1Ω2Mat)
            write(file, "AEmat",tabAEMat)
            write(file, "jELABmat",tabJMat)
            write(file, "tabvminmax",tabvminmax)
            write(file, "omgminmax",ωminmax)
        end

    end

end
