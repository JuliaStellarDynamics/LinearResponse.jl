


function GoStep(omgval::ComplexF64,
                identity_matrix::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                FHT::FiniteHilbertTransform.AbstractFHT,
                tabaMcoef::Array{Float64,4},
                tabωminωmax::Matrix{Float64},
                params::LinearParameters;
                ξ::Float64=1.0,
                DSIZE::Float64=1.e-4)

    # fill the M matrix
    tabM!(omgval,MMat,tabaMcoef,tabωminωmax,FHT,params)

    # compute the determinant of I-M
    permitivity_determinant_val = permitivity_determinant(identity_matrix,MMat,ξ=ξ)
    # if this value is bad, should immediately break

    # now compute a Jacobian via finite differences
    riomgval = omgval + DSIZE
    upomgval = omgval + DSIZE*1im

    tabM!(riomgval,MMat,tabaMcoef,tabωminωmax,FHT,params)
    permitivity_determinant_valri = permitivity_determinant(identity_matrix,MMat,ξ=ξ)

    tabM!(upomgval,MMat,tabaMcoef,tabωminωmax,FHT,params)
    permitivity_determinant_valup = permitivity_determinant(identity_matrix,MMat,ξ=ξ)

    dXirir = real(permitivity_determinant_valri-permitivity_determinant_val)/DSIZE
    dXirii = imag(permitivity_determinant_valri-permitivity_determinant_val)/DSIZE
    dXiupr = real(permitivity_determinant_valup-permitivity_determinant_val)/DSIZE
    dXiupi = imag(permitivity_determinant_valup-permitivity_determinant_val)/DSIZE

    # jacobian = [dXirir dXiupr ; dXirii dXiupi]
    # curpos = [real(permitivity_determinant_val);imag(permitivity_determinant_val)]
    # increment = jacobian \ (-curpos)

    increment1, increment2 = OrbitalElements._inverse_2d(dXirir,dXiupr,dXirii,dXiupi,-real(permitivity_determinant_val),-imag(permitivity_determinant_val))

    # omgval += increment[1] + increment[2]*1im
    omgval += increment1 + increment2*1im
    if params.VERBOSE > 0
        println("permitivity_determinant=$permitivity_determinant_val")
        println("omega=$omgval")
    end

    return omgval,permitivity_determinant_val
end


function FindDeterminantZero(startingomg::ComplexF64,
                             identity_matrix::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                             FHT::FiniteHilbertTransform.AbstractFHT,
                             tabaMcoef::Array{Float64,4},
                             tabωminωmax::Matrix{Float64},
                             params::LinearParameters;
                             ξ::Float64=1.0,
                             TOL::Float64=1.e-12,
                             DSIZE::Float64=1.e-4,
                             NSEARCH::Int64=100)

    # initial values
    permitivity_determinant_val = 1.0
    omgval = startingomg
    stepnum = 0

    while abs(permitivity_determinant_val) > TOL

        omgval,permitivity_determinant_val = GoStep(omgval,identity_matrix,MMat,FHT,tabaMcoef,tabωminωmax,params,ξ=ξ,DSIZE=DSIZE)

        if abs(permitivity_determinant_val) > 1.0
            break
        end

        stepnum += 1

        if stepnum > NSEARCH
            break
        end
        
    end

    # what happens when this goes wrong?
    if abs(permitivity_determinant_val) < TOL
        return omgval,permitivity_determinant_val
    end

    return NaN,NaN
end


function FindPole(startingomg::ComplexF64,
                  FHT::FiniteHilbertTransform.AbstractFHT,
                  params::LinearParameters;
                  ξ::Float64=1.0,
                  TOL::Float64=1.e-12,
                  DSIZE::Float64=1.e-4)

    # Preparinng computations of the response matrices
    MMat, tabaMcoef, tabωminωmax = PrepareM(params)
    identity_matrix = make_identity_matrix(params.nradial)

    bestomg,detval = FindDeterminantZero(startingomg,identity_matrix,MMat,FHT,tabaMcoef,tabωminωmax,params,ξ=ξ,TOL=TOL,DSIZE=DSIZE)

    (params.VERBOSE >= 0) && println("Best O for n1max=$(params.n1max),nradial=$(params.nradial) == $bestomg")

    return bestomg,detval
end
