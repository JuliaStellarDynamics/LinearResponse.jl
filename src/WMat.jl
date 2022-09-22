

"""
    Wintegrand(u,a,e)

Integrand computation/update for FT of basis elements
"""
function Wintegrand(w::Float64,
                    a::Float64,e::Float64,L::Float64,
                    Ω1::Float64,Ω2::Float64,
                    lharmonic::Int64,
                    basis::AstroBasis.Basis_type,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function;
                    EDGE::Float64=0.01,
                    VERBOSE::Int64=0)

    
    # Current location of the radius, r=r(w)
    rval = OrbitalElements.ru(w,a,e)
    
    # Current value of the radial frequency integrand (almost dθ/dw)
    gval = OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,w,a,e,EDGE=EDGE)

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,lharmonic,rval)

    # the velocity for integration (dθ1dw, dθ2dw)
    return Ω1*gval, (Ω2 - L/(rval^(2)))*gval 
end

"""
    WBasisFT(a,e,Ω1,Ω2,n1,n2,)

Fourier Transform of basis elements using RK4 scheme
result store in place
"""
function WBasisFT(a::Float64,e::Float64,
                  Ω1::Float64,Ω2::Float64,
                  n1::Int64,n2::Int64,
                  lharmonic::Int64,
                  basis::AstroBasis.Basis_type,
                  ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                  restab::Union{Array{Float64},SubArray{Float64}};
                  Kw::Int64=50,
                  EDGE::Float64=0.01,
                  VERBOSE::Int64=0)

    @assert length(restab) == basis.nmax "CallAResponse.WBasisFT: Result array not of the same size as the basis"

    # Integration step
    dw = (2.0)/(Kw)

    # need angular momentum
    Lval = OrbitalElements.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)

    # Initialise the state vectors: w, θ1, (θ2-psi)
    w, θ1, θ2 = -1.0, 0.0, 0.0

    # Initialize integrand
    dθ1dw, dθ2dw = Wintegrand(w,a,e,Lval,Ω1,Ω2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,EDGE=EDGE,VERBOSE=VERBOSE)

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for istep=1:Kw

        ####
        # RK4 Step 1
        ####
        # compute the first prefactor
        pref1 = (1.0/6.0)*dw*(1.0/(pi))*dθ1dw*cos(n1*θ1 + n2*θ2)

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            @inbounds restab[np] += pref1*basis.tabUl[np]
        end

        # update velocities at end of step 1
        dθ1_1 = dw*dθ1dw
        dθ2_1 = dw*dθ2dw

        ####
        # RK4 Step 2
        ####
        # Update the time by half a timestep
        w += 0.5*dw

        # Update integrand
        dθ1dw, dθ2dw = Wintegrand(w,a,e,Lval,Ω1,Ω2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,EDGE=EDGE,VERBOSE=VERBOSE)

        # common prefactor for all the increments
        # depends on the updated (θ1+0.5*dθ1_1,θ2+0.5*dθ2_1)
        # the factor (1.0/3.0) comes from RK4
        pref2 = (1.0/3.0)*dw*(1.0/(pi))*dθ1dw*cos(n1*(θ1+0.5*dθ1_1) + n2*(θ2+0.5*dθ2_1))

        # update velocities at end of step 2
        dθ1_2 = dw*dθ1dw
        dθ2_2 = dw*dθ2dw

        ####
        # RK4 step 3
        ####
        # The time, u, is not updated for this step
        # For this step, no need to re-compute the basis elements, as r has not been updated
        # Common prefactor for all the increments
        # Depends on the updated (θ1+0.5*dθ1_2,θ2+0.5*dθ2_2)
        # the factor (1.0/3.0) comes from RK4
        pref3 = (1.0/3.0)*dw*(1.0/(pi))*dθ1dw*cos(n1*(θ1+0.5*dθ1_2) + n2*(θ2+0.5*dθ2_2))

        # Loop over the radial indices to sum basis contributions
        # Contribution of steps 2 and 3 together
        for np=1:basis.nmax
            @inbounds restab[np] += (pref2+pref3)*basis.tabUl[np]
        end

        dθ1_3 = dθ1_2 # Does not need to be updated
        dθ2_3 = dθ2_2 # Does not need to be updated

        ####
        # RK4 step 4
        ####
        w += 0.5*dw # Updating the time by half a timestep: we are now at the next u value
        
        # Update integrand
        dθ1dw, dθ2dw = Wintegrand(w,a,e,Lval,Ω1,Ω2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,EDGE=EDGE,VERBOSE=VERBOSE)

        # Common prefactor for all the increments
        # Depends on the updated (θ1+dθ1_3,θ2+dθ2_3)
        # The factor (1.0/6.0) comes from RK4
        pref4 = (1.0/6.0)*dw*(1.0/(pi))*dθ1dw*cos(n1*(θ1+dθ1_3) + n2*(θ2+dθ2_3))

        # Loop over the radial indices to sum basis contributions
        for np=1:basis.nmax
            @inbounds restab[np] += pref4*basis.tabUl[np]
        end

        dθ1_4 = dw*dθ1dw # Current velocity for θ1
        dθ2_4 = dw*dθ2dw # Current velocity for θ2

        # Update the positions using RK4-like sum (Simpson's 1/3 rule)
        θ1 += (dθ1_1 + 2.0*dθ1_2 + 2.0*dθ1_3 + dθ1_4)/(6.0)
        θ2 += (dθ2_1 + 2.0*dθ2_2 + 2.0*dθ2_3 + dθ2_4)/(6.0)

        # clean or check nans?

    end # RK4 integration
end

"""
without Ω1, Ω2
"""
function WBasisFT(a::Float64,e::Float64,
                  n1::Int64,n2::Int64,
                  lharmonic::Int64,
                  basis::AstroBasis.Basis_type,
                  ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                  restab::Union{Array{Float64},SubArray{Float64}};
                  Kw::Int64=50,
                  EDGE::Float64=0.01,
                  TOLECC::Float64=0.001,
                  NINT::Int64=32,
                  VERBOSE::Int64=0)

    # Frequencies
    Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e;action=action,TOLECC=TOLECC,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)
    # Basis FT
    WBasisFT(a,e,Ω1,Ω2,n1,n2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,restab,Ω₀=Ω₀,Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)
end


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

            # Compute W(u,v) for every basis element using RK4 scheme
            WBasisFT(a,e,Ω1,Ω2,n1,n2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,view(tabWMat,:,kuval,kvval),Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)
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
    CheckConfigurationDirectories(wmatdir=wmatdir) || (return 0)

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
    (VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: Number of resonances to compute: $nbResVec")

    Threads.@threads for i = 1:nbResVec
        k = Threads.threadid()
        n1,n2 = tabResVec[1,i],tabResVec[2,i]

        (VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: Computing W for the ($n1,$n2) resonance.")

        # If it has been already computed
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
