
########################################################################
#
# Data structure
#
########################################################################

"""
@TO DESCRIBE
"""
function WMatdataCreateIsochrone(n1::Int64,n2::Int64,
                                 params::LinearParameters)

    # compute the frequency scaling factors for this resonance
    ωmin, ωmax = OrbitalElements.FindWminWmaxIsochrone(n1,n2)

    # Useful parameters
    nradial = params.nradial
    Ku, Kv = params.Ku, params.Kv

    return WMatdataType(ωmin,ωmax,zeros(Float64,2,Ku),
                        zeros(Float64,nradial,Kv,Ku),
                        zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku), # Orbital mappings
                        zeros(Float64,Kv,Ku))
end



########################################################################
#
# Fourier Transform at one location
#
########################################################################

"""
    WintegrandIsochrone(u,a,e)

Integrand computation/update for FT of basis elements
"""
function WintegrandIsochrone(w::Float64,
                             a::Float64,e::Float64,L::Float64,
                             Ω1::Float64,Ω2::Float64,
                             bc::Float64,M::Float64,G::Float64,
                             basis::AstroBasis.AbstractAstroBasis,
                             params::LinearParameters)::Tuple{Float64,Float64}

    # Current location of the radius, r=r(w)
    rval = OrbitalElements.ru(w,a,e)

    # Current value of the radial frequency integrand (almost dθ/dw)
    rp,ra = OrbitalElements.RpRaFromAE(a,e)
    gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=params.Orbitalparams.Ω₀)

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,params.lharmonic,rval)

    # the velocity for integration (dθ1dw, dθ2dw)
    return Ω1*gval, (Ω2 - L/(rval^(2)))*gval
end

"""
    WBasisFTIsochrone(a,e,Ω₁,Ω₂,n1,n2,bc,M,G,basisFT,restab,params)

Fourier Transform of basis elements using RK4 scheme
result stored in place
"""
function WBasisFTIsochrone(a::Float64,e::Float64,
                           Ω1::Float64,Ω2::Float64,
                           n1::Int64,n2::Int64,
                           bc::Float64,M::Float64,G::Float64,
                           basis::AstroBasis.AbstractAstroBasis,
                           restab::Vector{Float64},
                           params::LinearParameters)

    @assert length(restab) == basis.nradial "LinearResponse.WBasisFT: FT array not of the same size as the basis"

    # Integration step
    Kwp = (params.ADAPTIVEKW) ? ceil(Int64,params.Kw/(0.1+(1-e))) : params.Kw

    # Caution : Reverse integration (lower error at apocenter than pericenter)
    # -> Result to multiply by -1
    dw = (2.0)/(Kwp)

    # get (rp,ra)
    rp,ra = OrbitalElements.RpRaFromAE(a,e)

    # need angular momentum
    Lval = OrbitalElements.isochroneLfromrpra(rp,ra,bc,M,G)

    # Initialise the state vectors: u, theta1, (theta2-psi)
    u, theta1, theta2 = -1.0, 0.0, 0.0

    # launch the integration from the left boundary
    gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=params.Orbitalparams.Ω₀)

    # Current location of the radius, r=r(u)
    rval = OrbitalElements.ru(u,a,e)

    # the velocity for integration
    dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,params.lharmonic,rval)

    fill!(restab,0.0)

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for istep=1:Kwp

        # RK4 step 1
        # compute the first prefactor
        pref1 = (1.0/6.0)*dw*(1.0/(pi))*dt1du*cos(n1*theta1 + n2*theta2)

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nradial
            @inbounds restab[np] += pref1*basis.tabUl[np]
        end

        # update velocities at end of step 1
        k11 = dw*dt1du
        k21 = dw*dt2du

        # RK4 step 2
        # Update the time by half a timestep
        u += 0.5*dw

        # Current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=params.Orbitalparams.Ω₀)

        # Current value of dtheta1/du and dtheta2/du
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # recompute the basis functions for the changed radius value
        AstroBasis.tabUl!(basis,params.lharmonic,rval)

        # Common prefactor for all the increments
        # Depends on the updated (theta1+0.5*k11,theta2+0.5*k21)
        # the factor (1.0/3.0) comes from RK4
        pref2 = (1.0/3.0)*dw*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k11) + n2*(theta2+0.5*k21))

        # Loop over the radial indices to sum basis contributions
        #for np=1:basis.nradial
        #    tabWMat[np,kuval,kvval] += pref2*basis.tabUl[np]
        #end

        # update velocities at end of step 2
        k12 = dw*dt1du
        k22 = dw*dt2du

        # RK4 step 3
        # The time, u, is not updated for this step
        # For this step, no need to re-compute the basis elements, as r has not been updated

        # Common prefactor for all the increments
        # depends on the updated (theta1+0.5*k12,theta2+0.5*k22)
        # the factor (1.0/3.0) comes from RK4
        pref3 = (1.0/3.0)*dw*(1.0/(pi))*dt1du*cos(n1*(theta1+0.5*k12) + n2*(theta2+0.5*k22))


        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nradial
            #tabWMat[np,kuval,kvval] += pref3*basis.tabUl[np]
            @inbounds restab[np] += (pref2+pref3)*basis.tabUl[np]
        end

        k13 = k12 # Does not need to be updated
        k23 = k22 # Does not need to be updated

        # RK4 step 4
        # update the time by half a timestep: we are now at the next u value
        u += 0.5*dw

        # current location of the radius, r=r(u)
        rval = OrbitalElements.ru(u,a,e)

        # current value of Theta
        gval = OrbitalElements.ThetaRpRaIsochrone(rp,ra,u,bc=bc,Ω₀=params.Orbitalparams.Ω₀)

        # Current value of dtheta1/du and dtheta2/du, always well-posed
        dt1du, dt2du = Ω1*gval, (Ω2 - Lval/(rval^(2)))*gval

        # updated basis elements for new rval
        AstroBasis.tabUl!(basis,params.lharmonic,rval)

        # Common prefactor for all the increments:
        # Depends on the updated (theta1+k13,theta2+k23)
        # the factor (1.0/6.0) comes from RK4
        pref4 = (1.0/6.0)*dw*(1.0/(pi))*dt1du*cos(n1*(theta1+k13) + n2*(theta2+k23))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nradial
            @inbounds restab[np] += pref4*basis.tabUl[np]
        end

        # current velocities for theta1,theta2
        k14 = dw*dt1du
        k24 = dw*dt2du

        # update the positions using RK4-like sum
        theta1 += (k11 + 2.0*k12 + 2.0*k13 + k14)/(6.0)
        theta2 += (k21 + 2.0*k22 + 2.0*k23 + k24)/(6.0)

    end # RK4 integration

    return nothing
end


"""
with basisFT struct
"""
function WBasisFTIsochrone(a::Float64,e::Float64,
                  Ω1::Float64,Ω2::Float64,
                  n1::Int64,n2::Int64,
                  bc::Float64,M::Float64,G::Float64,
                  basisFT::BasisFTtype,
                  params::LinearParameters)

    # Basis FT
    WBasisFTIsochrone(a,e,Ω1,Ω2,n1,n2,bc,M,G,basisFT.basis,basisFT.UFT,params)
end

"""
without Ω1, Ω2
"""
function WBasisFTIsochrone(a::Float64,e::Float64,
                  n1::Int64,n2::Int64,
                  bc::Float64,M::Float64,G::Float64,
                  basisFT::BasisFTtype,
                  params::LinearParameters)

    # Frequencies
    Ω1, Ω2 = IsochroneOmega12FromAE(a,e,bc,M,G)

    # Basis FT
    WBasisFTIsochrone(a,e,Ω1,Ω2,n1,n2,bc,M,G,basisFT.basis,basisFT.UFT,params)
end


"""
    MakeWmatUVIsochrone(n1::Int64, n2::Int64, tabu::Vector{Float64}, basisFT::BasisFTtype, bc::Float64, M::Float64, G::Float64, params::LinearParameters)

Construct the matrix `Wdata` based on the given parameters for an isochrone potential.

# Arguments
- `n1::Int64`: Integer representing some parameter.
- `n2::Int64`: Integer representing some parameter.
- `tabu::Vector{Float64}`: Vector of `Float64` representing values of `u`.
- `basisFT::BasisFTtype`: Type representing some basis Fourier transform.
- `bc::Float64`: Float representing some parameter.
- `M::Float64`: Float representing mass.
- `G::Float64`: Float representing gravitational constant.
- `params::LinearParameters`: Type representing some linear parameters.

# Returns
- `Wdata`: A structure containing computed values.

# Description
`MakeWmatUVIsochrone` constructs the matrix `Wdata` based on the provided parameters. It iterates over each combination of `u` and `v` values, computes various orbital elements such as eccentricity, semi-major axis, frequencies, and more. It then computes the Jacobian of the transformation from `(α, β)` to `(E, L)` and saves it in `Wdata`. Finally, it computes `W(u, v)` for each basis element using a Runge-Kutta 4th order scheme and stores the result in `Wdata.tabW`.

"""
function MakeWmatUVIsochrone(n1::Int64,n2::Int64,
                             tabu::Vector{Float64},
                             basisFT::BasisFTtype,
                             bc::Float64,M::Float64,G::Float64,
                             params::LinearParameters)

    @assert length(tabu) == params.Ku "LinearResponse.WMat.MakeWmatUV: tabu length is not Ku."

    Orbitalparams = params.Orbitalparams
    Ω₀ = Orbitalparams.Ω₀

    # allocate the results matrices
    Wdata = WMatdataCreateIsochrone(n1,n2,params)
    ωmin, ωmax = Wdata.ωmin, Wdata.ωmax

    # start the loop
    for kuval = 1:params.Ku

        # get the current u value
        uval = tabu[kuval]

        (params.VERBOSE > 2) && println("LinearResponse.WMat.MakeWMat: on step $kuval of $Ku: u=$uval.")

        # get the corresponding v boundary values
        vmin,vmax = OrbitalElements.FindVminVmaxIsochrone(n1,n2,uval)

        # saving them
        Wdata.tabvminmax[1,kuval], Wdata.tabvminmax[2,kuval] = vmin, vmax

        # determine the step size in v
        δvp = 1.0/params.Kv


        for kvval = 1:params.Kv

            # get the current v value
            vp   = δvp*(kvval-0.5)
            vval = vFromvp(vp,vmin,vmax,params.VMAPN)

            ####
            # (u,v) -> (a,e)
            ####
            # (u,v) -> (α,β)
            α,β = OrbitalElements.αβFromUV(uval,vval,n1,n2,ωmin,ωmax)
            # (α,β) -> (Ω1,Ω2)
            Ω₁, Ω₂ = OrbitalElements.FrequenciesFromαβ(α,β,Ω₀)
            # (Ω1,Ω2) -> (a,e)
            #a,e = OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,Ω₁,Ω₂,Orbitalparams)
            a,e = OrbitalElements.IsochroneAEFromOmega1Omega2(Ω₁,Ω₂,bc,M,G)

            # get (rp,ra)
            rp,ra = OrbitalElements.RpRaFromAE(a,e)

            # need angular momentum
            Lval = OrbitalElements.isochroneLfromrpra(rp,ra,bc,M,G)

            (params.VERBOSE > 2) && print("v=$kvval,o1=$Ω1,o2=$Ω2;")

            # save (u,v) values for later
            Wdata.tabUV[1,kvval,kuval], Wdata.tabUV[2,kvval,kuval] = uval, vval
            # save (Ω1,Ω2) values for later
            Wdata.tabΩ1Ω2[1,kvval,kuval], Wdata.tabΩ1Ω2[2,kvval,kuval] = Ω₁, Ω₂
            # save (a,e) values for later
            Wdata.tabAE[1,kvval,kuval], Wdata.tabAE[2,kvval,kuval] = a, e
            # save (E,L) values for later
            Wdata.tabEL[1,kvval,kuval], Wdata.tabEL[2,kvval,kuval] = OrbitalElements.isochroneELfromαβ(α,β,bc,M,G)

            # compute the Jacobian of the (α,β) ↦ (E,L) mapping here. a little more expensive, but savings in the long run
            Wdata.tabJ[kvval,kuval] = OrbitalElements.IsochroneJacELtoαβ(α,β,bc,M,G)

            # Compute W(u,v) for every basis element using RK4 scheme
            WBasisFTIsochrone(a,e,Ω₁,Ω₂,n1,n2,bc,M,G,basisFT,params)

            for np = 1:basisFT.basis.nradial
                Wdata.tabW[np,kvval,kuval] = basisFT.UFT[np]
            end
        end
    end

    return Wdata
end


########################################################################
#
# Store FT on all (u,v) points for all resonances
#
########################################################################

"""
    RunWmat(ψ,dψ,d2ψ,FHT,basis[,params])

@TO DESCRIBE
"""
function RunWmatIsochrone(FHT::FiniteHilbertTransform.AbstractFHT,
                          bc::Float64,M::Float64,G::Float64,
                          basis::AstroBasis.AbstractAstroBasis,
                          params::LinearParameters)

    # check the directories + basis and FHT values against the Parameters
    CheckDirectories(params.wmatdir)
    CheckBasisCompatibility(basis,params)
    CheckFHTCompatibility(FHT,params)

    # FT bases prep.
    basisFT = BasisFTcreate(basis)
    basesFT = [deepcopy(basisFT) for k=1:Threads.nthreads()]

    # print the length of the list of resonance vectors
    (params.VERBOSE > 0) && println("LinearResponse.WMat.RunWmatIsochrone: Number of resonances to compute: $(params.nbResVec)")

    Threads.@threads for i = 1:params.nbResVec
        k = Threads.threadid()
        n1,n2 = params.tabResVec[1,i],params.tabResVec[2,i]

        (params.VERBOSE > 0) && println("LinearResponse.WMat.RunWmatIsochrone: Computing W for the ($n1,$n2) resonance.")

        # Output file name
        outputfilename = WMatFilename(n1,n2,params)

        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,params,"LinearResponse.WMat.RunWMat: ($n1,$n2) resonance") || continue

        # compute the W matrices in UV space: timing optional
        if (params.VERBOSE > 1) && (k == 1)
            @time Wdata = MakeWmatUVIsochrone(n1,n2,FHT.tabu,
                                     basesFT[k],bc,M,G,params)
        else
            Wdata = MakeWmatUVIsochrone(n1,n2,FHT.tabu,
                               basesFT[k],bc,M,G,params)
        end

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(outputfilename, "w") do file
            # Mappings parameters
            write(file, "omgmin",Wdata.ωmin)
            write(file, "omgmax",Wdata.ωmax)
            write(file, "tabvminmax",Wdata.tabvminmax)
            # Mappings
            write(file, "UVmat",Wdata.tabUV)
            write(file, "Omgmat",Wdata.tabΩ1Ω2)
            write(file, "AEmat",Wdata.tabAE)
            write(file, "ELmat",Wdata.tabEL)
            # Jacobians
            write(file, "jELABmat",Wdata.tabJ)
            # Basis FT
            write(file, "wmat",Wdata.tabW)
            # Parameters
            WriteParameters(file,params)
        end
    end
end
