


function GoStep(omgval::ComplexF64,
                IMat::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                FHT::FiniteHilbertTransform.AbstractFHT,
                tabaMcoef::Array{Float64,4},
                tabωminωmax::Matrix{Float64},
                params::LinearParameters;
                ξ::Float64=1.0,
                DSIZE::Float64=1.e-4)

    # fill the M matrix
    tabM!(omgval,MMat,tabaMcoef,tabωminωmax,FHT,params)

    # compute the determinant of I-M
    detXival = detXi(IMat,MMat,ξ=ξ)
    # if this value is bad, should immediately break

    # now compute a Jacobian via finite differences
    riomgval = omgval + DSIZE
    upomgval = omgval + DSIZE*1im

    tabM!(riomgval,MMat,tabaMcoef,tabωminωmax,FHT,params)
    detXivalri = detXi(IMat,MMat,ξ=ξ)

    tabM!(upomgval,MMat,tabaMcoef,tabωminωmax,FHT,params)
    detXivalup = detXi(IMat,MMat,ξ=ξ)

    dXirir = real(detXivalri-detXival)/DSIZE
    dXirii = imag(detXivalri-detXival)/DSIZE
    dXiupr = real(detXivalup-detXival)/DSIZE
    dXiupi = imag(detXivalup-detXival)/DSIZE

    # jacobian = [dXirir dXiupr ; dXirii dXiupi]
    # curpos = [real(detXival);imag(detXival)]
    # increment = jacobian \ (-curpos)

    increment1, increment2 = OrbitalElements.inverse2Dlinear(dXirir,dXiupr,dXirii,dXiupi,-real(detXival),-imag(detXival))

    # omgval += increment[1] + increment[2]*1im
    omgval += increment1 + increment2*1im
    if params.VERBOSE > 0
        println("detXi=$detXival")
        println("omega=$omgval")
    end

    return omgval,detXival
end


function FindDeterminantZero(startingomg::ComplexF64,
                             IMat::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                             FHT::FiniteHilbertTransform.AbstractFHT,
                             tabaMcoef::Array{Float64,4},
                             tabωminωmax::Matrix{Float64},
                             params::LinearParameters;
                             ξ::Float64=1.0,
                             TOL::Float64=1.e-12,
                             DSIZE::Float64=1.e-4,
                             NSEARCH::Int64=100)

    # initial values
    detXival = 1.0
    omgval = startingomg
    stepnum = 0

    while abs(detXival) > TOL

        omgval,detXival = GoStep(omgval,IMat,MMat,FHT,tabaMcoef,tabωminωmax,params,ξ=ξ,DSIZE=DSIZE)

        if abs(detXival) > 1.0
            break
        end

        stepnum += 1

        if stepnum > NSEARCH
            break
        end
        
    end

    # what happens when this goes wrong?
    if abs(detXival) < TOL
        return omgval,detXival
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
    IMat = makeIMat(params.nradial)

    bestomg,detval = FindDeterminantZero(startingomg,IMat,MMat,FHT,tabaMcoef,tabωminωmax,params,ξ=ξ,TOL=TOL,DSIZE=DSIZE)

    (params.VERBOSE >= 0) && println("Best O for n1max=$(params.n1max),nradial=$(params.nradial) == $bestomg")

    return bestomg,detval
end
