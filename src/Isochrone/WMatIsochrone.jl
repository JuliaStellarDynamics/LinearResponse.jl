"""The same as WMat.jl, but Isochrone ψ specific for testing

What's different:
-Use analytic inversion (IsochroneAEFromΩ1Ω2)
-Use analytic relationship between alpha and beta for circular orbits

"""

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform


"""MakeWmatUVIsochrone(n1,n2,Ku,Kv,lharmonic,basis[,Ω₀,Kw])

@IMPROVE: give basis a type?
@IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
@IMPROVE: parallelise by launching from both -1 and 1?
@IMPROVE: adaptively check for NaN values?

@WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl
"""
function MakeWmatUVIsochrone(n1::Int64,n2::Int64,
                             Kuvals::Matrix{Float64},
                             Kv::Int64,
                             lharmonic::Int64,
                             basis::AstroBasis.Basis_type,
                             Ω₀::Float64=1.,
                             Kw::Int64=20;
                             bc::Float64=1.0,M::Float64=1.0,G::Float64=1.0,EDGE::Float64=0.01)

    Ku = length(Kuvals)

    # compute the frequency scaling factors for this resonance
    # using exact values for the isochrone potential
    wmin,wmax = OrbitalElements.FindWminWmaxIsochrone(n1,n2)

    # allocate the results matrix
    tabWMat = zeros(basis.nmax,Ku,Kv)
    tabaMat = zeros(Ku,Kv)
    tabeMat = zeros(Ku,Kv)

    # set the matrix step size
    duWMat = (2.0)/(Kw)

    # start the loop
    for kuval in 1:Ku

        # get the current u value
        uval = Kuvals[kuval]

        # get the corresponding v values
        # exact values for the isochrone
        vmin,vmax = OrbitalElements.FindVminVmaxIsochrone(n1,n2,uval)

        # determine the step size in v
        deltav = (vmax - vmin)/(Kv)

        for kvval in 1:Kv

            # get the current v value
            vval = vmin + deltav*(kvval-0.5)

            # these are exact relationships
            alpha,beta = OrbitalElements.αβFromUV(uval,vval,n1,n2,wmin,wmax)
            omega1,omega2 = alpha*Ω₀,alpha*beta*Ω₀

            # big step: convert input (u,v) to (rp,ra)
            # convert from omega1,omega2 to (a,e) using exact Isochrone inversion
            a,e = OrbitalElements.IsochroneAEFromOmega1Omega2(omega1,omega2,bc,M,G)

            # save (a,e) values for later
            tabaMat[kuval,kvval] = a
            tabeMat[kuval,kvval] = e

            # get (rp,ra)
            rp,ra = OrbitalElements.RpRafromAE(a,e)

            # need angular momentum
            Lval = OrbitalElements.isochroneLfromrpra(rp,ra,bc,M,G)

            # Initialise the state vectors: u, theta1, (theta2-psi)
            u, theta1, theta2 = -1.0, 0.0, 0.0

            # launch the integration from the left boundary
            gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

            # Uses the Rozier 2019 notation for the mapping to u
            Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5

            # Current location of the radius, r=r(u): isn't this exactly rp?
            rval = OrbitalElements.ru(u,a,e)

            # the velocity for integration
            dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

            # collect the basis elements (in place!)
            AstroBasis.tabUl!(basis,lharmonic,rval)

            # start the integration loop now that we are initialised
            # at each step, we are performing an RK4-like calculation
            for istep=1:Kw

                # RK4 step 1
                # compute the first prefactor
                pref1 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref1*basis.tabUl[np]
                end

                # update velocities at end of step 1
                k11 = duWMat*dt1du
                k21 = duWMat*dt2du

                # RK4 step 2
                # Update the time by half a timestep
                u += 0.5*duWMat

                # Current location of the radius, r=r(u)
                rval = OrbitalElements.ru(u,a,e)

                # current value of Theta
                gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

                # Current value of dtheta1/du and dtheta2/du
                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

                # recompute the basis functions for the changed radius value
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # Common prefactor for all the increments
                # Depends on the updated (theta1+0.5*k11,theta2+0.5*k21)
                # the factor (1.0/3.0) comes from RK4
                pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k11) + n2*(theta2+0.5*k21))

                # Loop over the radial indices to sum basis contributions
                #for np=1:basis.nmax
                #    tabWMat[np,kuval,kvval] += pref2*basis.tabUl[np]
                #end

                # update velocities at end of step 2
                k12 = duWMat*dt1du
                k22 = duWMat*dt2du

                # RK4 step 3
                # The time, u, is not updated for this step
                # For this step, no need to re-compute the basis elements, as r has not been updated

                # Common prefactor for all the increments
                # depends on the updated (theta1+0.5*k12,theta2+0.5*k22)
                # the factor (1.0/3.0) comes from RK4
                pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k12) + n2*(theta2+0.5*k22))


                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    #tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
                    tabWMat[np,kuval,kvval] += (pref2+pref3)*basis.tabUl[np]
                end

                k13 = k12 # Does not need to be updated
                k23 = k22 # Does not need to be updated

                # RK4 step 4
                # update the time by half a timestep: we are now at the next u value
                u += 0.5*duWMat

                # current location of the radius, r=r(u)
                rval = OrbitalElements.ru(u,a,e)

                # current value of Theta
                gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

                # Current value of dtheta1/du and dtheta2/du, always well-posed
                dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval

                # updated basis elements for new rval
                AstroBasis.tabUl!(basis,lharmonic,rval)

                # Common prefactor for all the increments:
                # Depends on the updated (theta1+k13,theta2+k23)
                # the factor (1.0/6.0) comes from RK4
                pref4 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+k13) + n2*(theta2+k23))

                # Loop over the radial indices to sum basis contributions
                for np=1:basis.nmax
                    tabWMat[np,kuval,kvval] += pref4*basis.tabUl[np]
                end

                # current velocities for theta1,theta2
                k14 = duWMat*dt1du
                k24 = duWMat*dt2du

                # update the positions using RK4-like sum
                theta1 += (k11 + 2.0*k12 + 2.0*k13 + k14)/(6.0)
                theta2 += (k21 + 2.0*k22 + 2.0*k23 + k24)/(6.0)

            end # RK4 integration
        end
    end
    return tabWMat,tabaMat,tabeMat
end



function IntegrateOverOrbit!(Wvals::Array{Float64},T1::Float64,T2::Float64,
                             a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             lharmonic::Int64,
                             n1::Int64,n2::Int64,
                             basis::AstroBasis.Basis_type,
                             Kw::Int64=100;
                             bc::Float64=1.,Ω₀::Float64=1.0,M::Float64=1.0,G::Float64=1.0)


     # set the matrix step size
     duWMat = (2.0)/(Kw)

    # get (rp,ra)
    rp,ra = OrbitalElements.RpRafromAE(a,e)

    # need angular momentum
    Lval = OrbitalElements.isochroneLfromrpra(rp,ra,bc,M,G)

    println("rp=$rp,ra=$ra,Lval=$Lval")

    # Initialise the state vectors: u, theta1, (theta2-psi)
    u, theta1, theta2 = -1.0, 0.0, 0.0

    # launch the integration from the left boundary
    gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

    # Current location of the radius, r=r(u): isn't this exactly rp?
    rval = OrbitalElements.ru(u,a,e)

    # the velocity for integration
    dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,lharmonic,rval)

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for istep=1:Kw

        # RK4 step 1
        # compute the first prefactor
        pref1 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            Wvals[np] += pref1*basis.tabUl[np]
        end

        # update velocities at end of step 1
        k11 = duWMat*dt1du
        k21 = duWMat*dt2du

        # RK4 step 2
        # Update the time by half a timestep
        u += 0.5*duWMat

        # Current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

        # Current value of dtheta1/du and dtheta2/du
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # recompute the basis functions for the changed radius value
        AstroBasis.tabUl!(basis,lharmonic,rval)

        # Common prefactor for all the increments
        # Depends on the updated (theta1+0.5*k11,theta2+0.5*k21)
        # the factor (1.0/3.0) comes from RK4
        pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k11) + n2*(theta2+0.5*k21))

        # update velocities at end of step 2
        k12 = duWMat*dt1du
        k22 = duWMat*dt2du

        # RK4 step 3
        # The time, u, is not updated for this step
        # For this step, no need to re-compute the basis elements, as r has not been updated

        # Common prefactor for all the increments
        # depends on the updated (theta1+0.5*k12,theta2+0.5*k22)
        # the factor (1.0/3.0) comes from RK4
        pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k12) + n2*(theta2+0.5*k22))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            #tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
            Wvals[np] += (pref2+pref3)*basis.tabUl[np]
        end

        k13 = k12 # Does not need to be updated
        k23 = k22 # Does not need to be updated

        # RK4 step 4
        # update the time by half a timestep: we are now at the next u value
        u += 0.5*duWMat

        # current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

        # Current value of dtheta1/du and dtheta2/du, always well-posed
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # updated basis elements for new rval
        AstroBasis.tabUl!(basis,lharmonic,rval)

        # Common prefactor for all the increments:
        # Depends on the updated (theta1+k13,theta2+k23)
        # the factor (1.0/6.0) comes from RK4
        pref4 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+k13) + n2*(theta2+k23))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            Wvals[np] += pref4*basis.tabUl[np]
        end

        # current velocities for theta1,theta2
        k14 = duWMat*dt1du
        k24 = duWMat*dt2du

        # update the positions using RK4-like sum
        theta1 += (k11 + 2.0*k12 + 2.0*k13 + k14)/(6.0)
        theta2 += (k21 + 2.0*k22 + 2.0*k23 + k24)/(6.0)



    end

    T1 = theta1
    T2 = theta2

end




function IntegrateOverOrbitAllSteps!(Wvals::Array{Float64,2},
                                     uarr::Array{Float64},theta1arr::Array{Float64},theta2arr::Array{Float64},
                                     a::Float64,e::Float64,
                                     Ω1::Float64,Ω2::Float64,
                                     lharmonic::Int64,
                                     n1::Int64,n2::Int64,
                                     basis::AstroBasis.Basis_type,
                                     Kw::Int64=100;
                                     bc::Float64=1.,Ω₀::Float64=1.0,M::Float64=1.0,G::Float64=1.0)


    # set the matrix step size
    duWMat = (2.0)/(Kw)

    # get (rp,ra)
    rp,ra = OrbitalElements.RpRafromAE(a,e)

    # need angular momentum
    Lval = OrbitalElements.isochroneLfromrpra(rp,ra,bc,M,G)

    # Initialise the state vectors: u, theta1, (theta2-psi)
    u, theta1, theta2 = -1.0, 0.0, 0.0

    # launch the integration from the left boundary
    gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

    # Current location of the radius, r=r(u): isn't this exactly rp?
    rval = OrbitalElements.ru(u,a,e)

    # the velocity for integration
    dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,lharmonic,rval)

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for wstep=1:Kw

        # RK4 step 1
        # compute the first prefactor
        pref1 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            Wvals[np,wstep] = pref1*basis.tabUl[np]
        end

        # update velocities at end of step 1
        k11 = duWMat*dt1du
        k21 = duWMat*dt2du

        # RK4 step 2
        # Update the time by half a timestep
        u += 0.5*duWMat

        # Current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

        # Current value of dtheta1/du and dtheta2/du
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # recompute the basis functions for the changed radius value
        AstroBasis.tabUl!(basis,lharmonic,rval)

        # Common prefactor for all the increments
        # Depends on the updated (theta1+0.5*k11,theta2+0.5*k21)
        # the factor (1.0/3.0) comes from RK4
        pref2 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k11) + n2*(theta2+0.5*k21))

        # update velocities at end of step 2
        k12 = duWMat*dt1du
        k22 = duWMat*dt2du

        # RK4 step 3
        # The time, u, is not updated for this step
        # For this step, no need to re-compute the basis elements, as r has not been updated

        # Common prefactor for all the increments
        # depends on the updated (theta1+0.5*k12,theta2+0.5*k22)
        # the factor (1.0/3.0) comes from RK4
        pref3 = (1.0/3.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k12) + n2*(theta2+0.5*k22))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            #tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
            Wvals[np,wstep] += (pref2+pref3)*basis.tabUl[np]
        end

        k13 = k12 # Does not need to be updated
        k23 = k22 # Does not need to be updated

        # RK4 step 4
        # update the time by half a timestep: we are now at the next u value
        u += 0.5*duWMat

        # current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=Ω₀)

        # Current value of dtheta1/du and dtheta2/du, always well-posed
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # updated basis elements for new rval
        AstroBasis.tabUl!(basis,lharmonic,rval)

        # Common prefactor for all the increments:
        # Depends on the updated (theta1+k13,theta2+k23)
        # the factor (1.0/6.0) comes from RK4
        pref4 = (1.0/6.0)*duWMat*(1.0/(pi))*dt1du*cos(n1*(theta1+k13) + n2*(theta2+k23))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            Wvals[np,wstep] += pref4*basis.tabUl[np]
        end

        # current velocities for theta1,theta2
        k14 = duWMat*dt1du
        k24 = duWMat*dt2du

        # update the positions using RK4-like sum
        theta1 += (k11 + 2.0*k12 + 2.0*k13 + k14)/(6.0)
        theta2 += (k21 + 2.0*k22 + 2.0*k23 + k24)/(6.0)

        # record values
        theta1arr[wstep] = (k11 + 2.0*k12 + 2.0*k13 + k14)/(6.0)
        theta2arr[wstep] = (k21 + 2.0*k22 + 2.0*k23 + k24)/(6.0)
        uarr[wstep] = u - 0.5*duWMat

    end

end




"""
    RunWmatIsochrone()

"""
function RunWmatIsochrone(wmatdir::String,
                          Ku::Int64,Kv::Int64,Kw::Int64,
                          basis::AstroBasis.Basis_type,
                          lharmonic::Int64,
                          n1max::Int64,
                          nradial::Int64,
                          Ω₀::Float64,
                          modelname::String,
                          rb::Float64;
                          bc::Float64=1.0,G::Float64=1.0,M::Float64=1.0,
                          VERBOSE::Int64=0,
                          OVERWRITE::Bool=false)

    # check wmat directory before proceeding (save time if not.)
    CheckConfigurationDirectories([wmatdir]) || (return 0)

    # get basis parameters
    ndim = basis.dimension
    nradialmax = basis.nmax

    # check if we can cover the specified radial orders
    if nradialmax > nradial
        println("CallAResponse.WMat.RunWmatIsochrone: the input basis does not have sufficient nradial ($nradialmax) for the requested value ($nradial).")
        return 0
    end

    # bases prep.
    bases=[deepcopy(basis) for k=1:Threads.nthreads()]

    # Legendre integration prep.
    tabuGLquadtmp,tabwGLquad = FiniteHilbertTransform.tabuwGLquad(Ku)
    tabuGLquad = reshape(tabuGLquadtmp,Ku,1)

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    # print the length of the list of resonance vectors
    println("CallAResponse.WMatIsochrone.RunWmatIsochrone: Number of resonances to compute: $nbResVec")

    Threads.@threads for i = 1:nbResVec
        k = Threads.threadid()
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        (VERBOSE > 0) && println("CallAResponse.WMatIsochrone.RunWmatIsochrone: Computing W for the ($n1,$n2) resonance.")

        if isfile(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw))
            file = h5open(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw), "r")
            oldnradial = read(file,"nradial")
            if (OVERWRITE == false) && (nradial <= oldnradial)
                (VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with higher nradial: no computation.")
                continue
            else
                (VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with lower nradial: recomputing and overwritting.")
            end
            close(file)
        end


        # currently defaulting to timed version:
        # could make this a flag (timing optional)
        if VERBOSE > 1
            @time tabWMat,tabaMat,tabeMat = MakeWmatUVIsochrone(n1,n2,tabuGLquad,Kv,lharmonic,bases[k],Ω₀,Kw,bc=bc,M=M,G=G)
        else
            tabWMat,tabaMat,tabeMat = MakeWmatUVIsochrone(n1,n2,tabuGLquad,Kv,lharmonic,bases[k],Ω₀,Kw,bc=bc,M=M,G=G)
        end

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(WMatFilename(wmatdir,modelname,lharmonic,n1,n2,rb,Ku,Kv,Kw), "w") do file
            write(file, "nradial",nradial)
            write(file, "wmat",tabWMat)
        end

    end

end
