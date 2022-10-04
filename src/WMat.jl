

"""
structure to store Fourier Transform values of basis elements
"""
struct BasisFT_type
    basis::AstroBasis.Basis_type
    UFT::Array{Float64}
end

function BasisFT_create(basis::AstroBasis.Basis_type)

    return BasisFT_type(basis,zeros(Float64,basis.nmax))
end

"""
    Wintegrand(u,a,e)

Integrand computation/update for FT of basis elements
"""
function Wintegrand(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    w::Float64,
                    a::Float64,e::Float64,L::Float64,
                    Ω1::Float64,Ω2::Float64,
                    basis::AstroBasis.Basis_type,
                    Parameters::ResponseParameters)


    # Current location of the radius, r=r(w)
    rval = OrbitalElements.ru(w,a,e)

    # Current value of the radial frequency integrand (almost dθ/dw)
    gval = OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,w,a,e,EDGE=Parameters.EDGE)

    # collect the basis elements (in place!)
    AstroBasis.tabUl!(basis,Parameters.lharmonic,rval)

    # the velocity for integration (dθ1dw, dθ2dw)
    return Ω1*gval, (Ω2 - L/(rval^(2)))*gval
end

"""
    WBasisFT(a,e,Ω1,Ω2,n1,n2,)

Fourier Transform of basis elements using RK4 scheme
result stored in place
"""
function WBasisFT(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                  a::Float64,e::Float64,
                  Ω1::Float64,Ω2::Float64,
                  n1::Int64,n2::Int64,
                  basis::AstroBasis.Basis_type,
                  Parameters::ResponseParameters)

    # Integration step
    dw = (2.0)/(Kw)

    # need angular momentum
    Lval = OrbitalElements.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)

    # Initialise the state vectors: w, θ1, (θ2-psi)
    w, θ1, θ2 = -1.0, 0.0, 0.0

    # Initialize integrand
    basis = basisFT.basis
    dθ1dw, dθ2dw = Wintegrand(ψ,dψ,d2ψ,d3ψ,w,a,e,Lval,Ω1,Ω2,Parameters)

    # Initialize container
    restab = basisFT.UFT
    fill!(restab,0.0)
    @assert length(restab) == basis.nmax "CallAResponse.WBasisFT: FT array not of the same size as the basis"

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for istep=1:Parameters.Kw

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
        dθ1dw, dθ2dw = Wintegrand(ψ,dψ,d2ψ,d3ψ,w,a,e,Lval,Ω1,Ω2,Parameters)

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
        dθ1dw, dθ2dw = Wintegrand(ψ,dψ,d2ψ,d3ψ,w,a,e,Lval,Ω1,Ω2,Parameters)

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
                  basisFT::BasisFT_type,
                  ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                  Parameters::ResponseParameters
                  )

    # Frequencies
    Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=Parameters.TOLECC,VERBOSE=Parameters.VERBOSE,NINT=Parameters.NINT,EDGE=Parameters.EDGE)
    # Basis FT
    WBasisFT(ψ,dψ,d2ψ,d3ψ,a,e,Ω1,Ω2,n1,n2,basisFT,Parameters)
end

"""
structure to store the W matrix computation results
"""
struct WMatdata_type
    tabW::Array{Float64,3}
    tabUV::Array{Float64,3}
    tabΩ1Ω2::Array{Float64,3}
    tabAE::Array{Float64,3}
    tabEL::Array{Float64,3}
    tabJ::Array{Float64,2}
    ωminmax::Vector{Float64}
    tabvminmax::Array{Float64,2}
end

"""
@TO DESCRIBE
"""
function WMatdata_create(nmax::Int64,Ku::Int64,Kv::Int64)

    return WMatdata_type(zeros(Float64,nmax,Ku,Kv),
                         zeros(Float64,Ku,Kv,2),zeros(Float64,Ku,Kv,2),zeros(Float64,Ku,Kv,2),zeros(Float64,Ku,Kv,2), # Orbital mappings
                         zeros(Float64,Ku,Kv),
                         zeros(Float64,2),zeros(Float64,Ku,2))
end

"""
@TO DESCRIBE: from known data
"""
function WMatdata_create(tabW::Array{Float64,3},
                         tabUV::Array{Float64,3},tabΩ1Ω2::Array{Float64,3},tabAE::Array{Float64,3},tabEL::Array{Float64,3},
                         tabJ::Array{Float64,2},
                         ωminmax::Vector{Float64},tabvminmax::Array{Float64,2})

    return WMatdata_type(tabW,
                         tabUV,tabΩ1Ω2,tabAE,tabEL,
                         tabJ,
                         ωminmax,tabvminmax)
end

"""
    MakeWmatUV(ψ,dψ,d2ψ,d3ψ,n1,n2,basisFT,Parameters)

@TO DESCRIBE
"""
function MakeWmatUV(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,βc::Function,
                    n1::Int64,n2::Int64,
                    tabu::Vector{Float64},
                    basisFT::BasisFT_type,
                    Parameters::ResponseParameters)

    """
    @IMPROVE: consolidate steps 2 and 3, which only have different prefactors from velocities
    @IMPROVE: parallelise by launching from both -1 and 1?
    @IMPROVE: adaptively check for NaN values?

    @WARNING: when parallelising, basis will need to be copies so it can overwrite tabUl
    """

    # get the number of u samples from the input vector of u vals
    Ku = length(tabu)

    # allocate the results matrices
    Wdata = WMatdata_create(basisFT.basis.nmax,Ku,Kv)

    # Frequency cuts associated to [rmin,rmax]
    # @IMPROVE: compute them once (independant of n1,n2) and function argument ?
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)
    # compute the frequency scaling factors for this resonance
    ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
    Wdata.ωminmax[1], Wdata.ωminmax[2] = ωmin, ωmax

    # start the loop
    for kuval in 1:Ku

        # get the current u value
        uval = tabu[kuval]

        (Parameters.VERBOSE > 2) && println("CallAResponse.WMat.MakeWMat: on step $kuval of $Ku: u=$uval.")

        # get the corresponding v boundary values
        vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

        # saving them
        Wdata.tabvminmax[kuval,1], Wdata.tabvminmax[kuval,2] = vmin, vmax

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
            a,e = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,NINT=Parameters.NINT,EDGE=Parameters.EDGE,VERBOSE=Parameters.VERBOSE)

            (Parameters.VERBOSE > 2) && print("v=$kvval,o1=$Ω1,o2=$Ω2;")

            # save (u,v) values for later
            Wdata.tabUV[kuval,kvval,1], Wdata.tabUV[kuval,kvval,2] = uval, vval
            # save (Ω1,Ω2) values for later
            Wdata.tabΩ1Ω2[kuval,kvval,1], Wdata.tabΩ1Ω2[kuval,kvval,2] = Ω1, Ω2
            # save (a,e) values for later
            Wdata.tabAE[kuval,kvval,1], Wdata.tabAE[kuval,kvval,2] = a, e

            # compute the Jacobian (E,L)->(alpha,beta) here. a little more expensive, but savings in the long run
            Wdata.tabJ[kuval,kvval] = OrbitalElements.JacELToαβAE(a,e,ψ,dψ,d2ψ,Parameters.Ω₀)

            # Compute W(u,v) for every basis element using RK4 scheme
            WBasisFT(ψ,dψ,d2ψ,d3ψ,a,e,Ω1,Ω2,n1,n2,basisFT,Parameters)
            for np = 1:basisFT.basis.nmax
                Wdata.tabW[np,kuval,kvval] = basisFT.UFT[np]
            end
        end
    end

    return Wdata
end


"""
    RunWmat(ψ,dψ,d2ψ,d3ψ,wmatdir,Ku,Kv,Kw,basis,lharmonic,n1max,nradial,Ω₀,modelname,rb,rmin,rmax[,VERBOSE,OVERWRITE])

@TO DESCRIBE
"""
function RunWmat(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                 FHT::FiniteHilbertTransform.FHTtype,
                 basis::AstroBasis.Basis_type,
                 Parameters::ResponseParameters)

    # check wmat directory before proceeding (save time if not.)
    CheckConfigurationDirectories([Parameters.wmatdir]) || (return 0)

    # get basis parameters
    ndim, nradial, rb = basis.dimension, basis.nmax, basis.rb

    # bases prep.
    AstroBasis.fill_prefactors!(basis)
    basisFT = BasisFT_create(basis)
    basesFT=[deepcopy(basisFT) for k=1:Threads.nthreads()]

    # define a function for βcircular
    βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin=Parameters.rmin,rmax=Parameters.rmax)

    # Integration points
    tabu, Ku = FHT.tabu, FHT.Ku

    # Resonance vectors
    #nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)
    #Parameters.nbResVec
    #Parameters.tabResVec

    # print the length of the list of resonance vectors
    (Parameters.VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: Number of resonances to compute: $(Parameters.nbResVec)")

    Threads.@threads for i = 1:Parameters.nbResVec
        k = Threads.threadid()
        n1,n2 = Parameters.tabResVec[1,i],Parameters.tabResVec[2,i]

        (Parameters.VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: Computing W for the ($n1,$n2) resonance.")

        # If it has been already computed
        if isfile(WMatFilename(n1,n2,Parameters))
            file = h5open(WMatFilename(n1,n2,Parameters), "r")
            oldnradial = read(file,"nradial")
            if (Parameters.OVERWRITE == false) && (Parameters.nradial <= oldnradial)
                (Parameters.VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with higher nradial: no computation.")
                continue
            else
                (Parameters.VERBOSE > 0) && println("CallAResponse.WMat.RunWmat: ($n1,$n2) resonanance WMat file already exists with lower nradial: recomputing and overwritting.")
            end
            close(file)
        end

        # compute the W matrices in UV space: timing optional
        if (Parameters.VERBOSE > 1)
            @time Wdata = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,βc,
                                     n1,n2,FHT.tabu,
                                     basesFT[k],Parameters)
        else
            Wdata = MakeWmatUV(ψ,dψ,d2ψ,d3ψ,βc,
                               n1,n2,FHT.tabu,
                               basesFT[k],Parameters)
        end

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(WMatFilename(n1,n2,Parameters), "w") do file
            write(file, "nradial",nradial)
            write(file, "wmat",Wdata.tabW)
            write(file, "Omgmat",Wdata.tabΩ1Ω2)
            write(file, "AEmat",Wdata.tabAE)
            write(file, "jELABmat",Wdata.tabJ)
            write(file, "tabvminmax",Wdata.tabvminmax)
            write(file, "omgminmax",Wdata.ωminmax)
        end
    end
end
