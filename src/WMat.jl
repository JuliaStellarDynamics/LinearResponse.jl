
########################################################################
#
# Fourier Transform over angles for a given orbit
#
########################################################################
"""
    angle_fouriertransform!(
        result::Vector{Float64},
        fun::Function,
        a::Float64,
        e::Float64,
        n1::Int64,
        n2::Int64,
        model::Potential,
        params::LinearParameters;
        L::Float64=0.0,
        Ω1::Float64=0.0,
        Ω2::Float64=0.0
    )

Perform the angle Fourier transform of a given function `fun` using the specified parameters.
The result is stored in the `result` vector.

Assuming the function is of the form `fun(r)`, the Fourier transform is given by:  
    `result = ∫[0,π] fun(r) * cos(n1*θ1 + n2*(θ2-ϕ)) dθ1`
as `r` is an even function of `θ1` alone.

# Arguments
- `result::Vector{Float64}`: The vector to store the result of the Fourier transform.
- `fun::Function`: The function to be transformed. Has to be a function of radius only.
- `a::Float64`: Semi-major axis.
- `e::Float64`: Eccentricity.
- `n1::Int64`: Fourier mode for the first angle.
- `n2::Int64`: Fourier mode for the second angle.
- `model::Potential`: The potential model.
- `params::LinearParameters`: Additional parameters for the transformation.

# Optional Arguments
- `L::Float64=0.0`: Angular momentum. If not provided, it will be calculated 
from `a`, `e`, `model`, and `params`.
- `Ω1::Float64=0.0`: Radial frequency. If not provided, it will be calculated 
from `a`, `e`, `model`, and `params`.
- `Ω2::Float64=0.0`: Azimuthal frequency. If not provided, it will be calculated 
from `a`, `e`, `model`, and `params`.
"""
function angle_fouriertransform!(
    result::Vector{Float64},
    fun::F0,
    a::Float64,
    e::Float64,
    n1::Int64,
    n2::Int64,
    model::Potential,
    params::LinearParameters;
    L::Float64=0.0,
    Ω1::Float64=0.0,
    Ω2::Float64=0.0
) where {F0<:Function}
    output_length = length(fun(0.0))
    @assert length(result) == output_length "Result array length does not 
        match the function's output length."

    # Integration step
    Kwp = (params.ADAPTIVEKW) ? ceil(Int64,params.Kw/(0.1+(1-e))) : params.Kw
    Orbitalparams = params.Orbitalparams
    
    if L == 0.0
        # need angular momentum
        _, L = EL_from_ae(a, e, model, Orbitalparams)
    end
    if Ω1 == 0.0 || Ω2 == 0.0
        # need frequencies
        Ω1, Ω2 = frequencies_from_ae(a, e, model, Orbitalparams)
    end

    # Caution : Reverse integration (lower error at apocenter than pericenter)
    # -> Result to multiply by -1
    dw = -2/Kwp

    # Initialise the state vectors: w, θ1, (θ2-psi)
    # @WARNING: in this function θ2 really stands for (θ2-psi) which is a function 
    # of the anomaly `w`.
    # Reverse integration, starting at apocenter
    w, θ1, θ2 = 1.0, pi, 0.0

    # Initialize integrand
    dθ1dw, dθ2dw = angles_gradient(w, a, e, model, Orbitalparams, L=L, Ω1=Ω1, Ω2=Ω2)
    integrand = fun(radius_from_anomaly(w, a, e, model, Orbitalparams))

    # Initialize container
    fill!(result,0.0)

    # start the integration loop now that we are initialised
    # at each step, we are performing an RK4-like calculation
    for step in 1:Kwp

        ####
        # RK4 Step 1
        ####
        # compute the first prefactor
        pref1 = (1/6) * dw * (1/pi) * dθ1dw * cos(n1*θ1 + n2*θ2)

        # Add contribution to the result
        result .+= pref1 .* integrand

        # update velocities at end of step 1
        dθ1_1 = dw * dθ1dw
        dθ2_1 = dw * dθ2dw

        ####
        # RK4 Step 2
        ####
        # Update the time by half a timestep
        w += 0.5*dw

        # Update integrand
        dθ1dw, dθ2dw = angles_gradient(w, a, e, model, Orbitalparams, L=L, Ω1=Ω1, Ω2=Ω2)
        integrand = fun(radius_from_anomaly(w, a, e, model, Orbitalparams))

        # common prefactor for all the increments
        # depends on the updated (θ1+0.5*dθ1_1,θ2+0.5*dθ2_1)
        # the factor 1/3 comes from RK4
        pref2 = (1/3) * dw * (1/pi) * dθ1dw * cos(n1*(θ1+0.5*dθ1_1) + n2*(θ2+0.5*dθ2_1))

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
        # the factor 1/3 comes from RK4
        pref3 = (1/3) * dw * (1/pi) * dθ1dw * cos(n1*(θ1+0.5*dθ1_2) + n2*(θ2+0.5*dθ2_2))

        # Add contribution of steps 2 and 3 together
        result .+= (pref2 + pref3) .* integrand

        dθ1_3 = dθ1_2 # Does not need to be updated
        dθ2_3 = dθ2_2 # Does not need to be updated

        ####
        # RK4 step 4
        ####
        w += 0.5*dw # Updating the time by half a timestep: we are now at the next u value

        # Update integrand
        dθ1dw, dθ2dw = angles_gradient(w, a, e, model, Orbitalparams, L=L, Ω1=Ω1, Ω2=Ω2)
        integrand = fun(radius_from_anomaly(w, a, e, model, Orbitalparams))

        # Common prefactor for all the increments
        # Depends on the updated (θ1+dθ1_3,θ2+dθ2_3)
        # The factor (1/6) comes from RK4
        pref4 = (1/6) * dw * (1/pi) * dθ1dw * cos(n1*(θ1+dθ1_3) + n2*(θ2+dθ2_3))

        # Loop over the radial indices to sum basis contributions
        result .+= pref4 .* integrand

        dθ1_4 = dw*dθ1dw # Current velocity for θ1
        dθ2_4 = dw*dθ2dw # Current velocity for θ2

        # Update the positions using RK4-like sum (Simpson's 1/3 rule)
        θ1 += (dθ1_1 + 2*dθ1_2 + 2*dθ1_3 + dθ1_4)/6
        θ2 += (dθ2_1 + 2*dθ2_2 + 2*dθ2_3 + dθ2_4)/6

        # clean or check nans?

    end # RK4 integration

    # -1 factor (reverse integration)
    result .*= -1

    return result
end

function angle_fouriertransform(
    fun::F0,
    a::Float64,
    e::Float64,
    n1::Int64,
    n2::Int64,
    model::Potential,
    params::LinearParameters;
    L::Float64=0.0,
    Ω1::Float64=0.0,
    Ω2::Float64=0.0
) where {F0<:Function}
    result = zeros(Float64,length(fun(0.0)))
    return angle_fouriertransform!(
        result, fun, a, e, n1, n2, model, params, L=L, Ω1=Ω1, Ω2=Ω2
    )
end

"""
for a basis element, compute the Fourier transform over the angles.
"""
function angle_fouriertransform(
    basis::AstroBasis.AbstractAstroBasis,
    a::Float64,
    e::Float64,
    n1::Int64,
    n2::Int64,
    model::Potential,
    params::LinearParameters;
    L::Float64=0.0,
    Ω1::Float64=0.0,
    Ω2::Float64=0.0
)
    function basisft_integrand(r::Float64)
        # collect the basis elements (in place!)
        AstroBasis.tabUl!(basis, params.lharmonic, r)
        return basis.tabUl
    end
    return angle_fouriertransform(
        basisft_integrand, a, e, n1, n2, model, params, L=L, Ω1=Ω1, Ω2=Ω2
    )
end

########################################################################
#
# Data structure
#
########################################################################
"""
structure to store the W matrix computation results
"""
struct FourierTransformedBasisData
    ωmin::Float64
    ωmax::Float64
    Ω₀::Float64
    tabvminmax::Array{Float64,2}

    tabW::Array{Float64,3}
    tabUV::Array{Float64,3}
    tabΩ1Ω2::Array{Float64,3}
    tabAE::Array{Float64,3}
    tabEL::Array{Float64,3}
    tabJ::Array{Float64,2}
end

"""
@TO DESCRIBE
"""
function WMatdataCreate(model::OrbitalElements.Potential,
                        n1::Int64,n2::Int64,
                        params::LinearParameters)

    # compute the frequency scaling factors for this resonance
    ωmin, ωmax = frequency_extrema(n1,n2,model,params.Orbitalparams)

    # Useful parameters
    nradial = params.nradial
    Ku, Kv = params.Ku, params.Kv

    return FourierTransformedBasisData(ωmin,ωmax,frequency_scale(model),zeros(Float64,2,Ku),
                        zeros(Float64,nradial,Kv,Ku),
                        zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku),zeros(Float64,2,Kv,Ku), # Orbital mappings
                        zeros(Float64,Kv,Ku))
end

########################################################################
#
# Store FT on all (u,v) points for one resonance number
#
########################################################################

"""
@TO DESCRIBE
"""
function MakeWmatUV(model::OrbitalElements.Potential,
                    n1::Int64,n2::Int64,
                    tabu::Vector{Float64},
                    basis::AstroBasis.AbstractAstroBasis,
                    params::LinearParameters)

    @assert length(tabu) == params.Ku "LinearResponse.WMat.MakeWmatUV: tabu length is not Ku."

    Orbitalparams = params.Orbitalparams
    Ω₀ = frequency_scale(model)

    # allocate the results matrices
    Wdata = WMatdataCreate(model,n1,n2,params)

    # start the loop
    for kuval = 1:params.Ku

        # get the current u value
        uval = tabu[kuval]

        (params.VERBOSE > 2) && println("LinearResponse.WMat.MakeWMat: on step $kuval of $Ku: u=$uval.")

        # get the corresponding v boundary values
        resonance = Resonance(n1,n2,model,Orbitalparams)
        vmin,vmax = v_boundaries(uval,resonance,model,Orbitalparams)

        # saving them
        Wdata.tabvminmax[1,kuval], Wdata.tabvminmax[2,kuval] = vmin, vmax

        # determine the step size in v
        δvp = 1.0/params.Kv


        for kvval = 1:params.Kv

            # get the current v value
            vp   = δvp*(kvval-0.5)
            vval = v_from_vp(vp, vmin, vmax, n=params.VMAPN)

            ####
            # (u,v) -> (a,e)
            ####
            # (u,v) -> (α,β)
            α,β = αβ_from_uv(uval,vval,resonance)
            # (α,β) -> (Ω1,Ω2)
            Ω₁, Ω₂ = frequencies_from_αβ(α,β,Ω₀)
            # (Ω1,Ω2) -> (a,e)
            a,e = ae_from_frequencies(Ω₁,Ω₂,model,Orbitalparams)

            (params.VERBOSE > 2) && print("v=$kvval,o1=$Ω1,o2=$Ω2;")

            # save (u,v) values for later
            Wdata.tabUV[1,kvval,kuval], Wdata.tabUV[2,kvval,kuval] = uval, vval
            # save (Ω1,Ω2) values for later
            Wdata.tabΩ1Ω2[1,kvval,kuval], Wdata.tabΩ1Ω2[2,kvval,kuval] = Ω₁, Ω₂
            # save (a,e) values for later
            Wdata.tabAE[1,kvval,kuval], Wdata.tabAE[2,kvval,kuval] = a, e
            # save (E,L) values for later
            E, L = EL_from_ae(a,e,model,Orbitalparams)
            Wdata.tabEL[1,kvval,kuval], Wdata.tabEL[2,kvval,kuval] = E, L

            # compute the Jacobian of the (α,β) ↦ (E,L) mapping here. a little more expensive, but savings in the long run
            Wdata.tabJ[kvval,kuval] = ae_to_EL_jacobian(a,e,model,Orbitalparams)/ae_to_αβ_jacobian(a,e,model,Orbitalparams)

            # Compute W(u,v) for every basis element using RK4 scheme
            basisft = angle_fouriertransform(
                basis, a, e, n1, n2, model, params, L=L, Ω1=Ω₁, Ω2=Ω₂
            )

            for np = 1:basis.nradial
                Wdata.tabW[np,kvval,kuval] = basisft[np]
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
function RunWmat(model::OrbitalElements.Potential,
                 FHT::FiniteHilbertTransform.AbstractFHT,
                 basis::AstroBasis.AbstractAstroBasis,
                 params::LinearParameters)

    # check the directories + basis and FHT values against the Parameters
    CheckDirectories(params.wmatdir)
    CheckBasisCompatibility(basis,params)
    CheckFHTCompatibility(FHT,params)

    # FT bases prep.
    bases = [deepcopy(basis) for k=1:Threads.nthreads()]

    # print the length of the list of resonance vectors
    (params.VERBOSE > 0) && println("LinearResponse.WMat.RunWmat: Number of resonances to compute: $(params.nbResVec)")

    Threads.@threads for i = 1:params.nbResVec
        k = Threads.threadid()
        n1,n2 = params.tabResVec[1,i],params.tabResVec[2,i]

        (params.VERBOSE > 0) && println("LinearResponse.WMat.RunWmat: Computing W for the ($n1,$n2) resonance.")

        # Output file name
        outputfilename = WMatFilename(n1,n2,params)

        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,params,"LinearResponse.WMat.RunWMat: ($n1,$n2) resonance") || continue

        # compute the W matrices in UV space: timing optional
        if (params.VERBOSE > 1) && (k == 1)
            @time Wdata = MakeWmatUV(model,
                                     n1,n2,FHT.tabu,
                                     bases[k],params)
        else
            Wdata = MakeWmatUV(model,
                               n1,n2,FHT.tabu,
                               bases[k],params)
        end

        # now save: we are saving not only W(u,v), but also a(u,v) and e(u,v).
        # could consider saving other quantities as well to check mappings.
        h5open(outputfilename, "w") do file
            # Mappings parameters
            write(file, "omgmin",Wdata.ωmin)
            write(file, "omgmax",Wdata.ωmax)
            write(file,"Ω₀",Wdata.Ω₀)
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
