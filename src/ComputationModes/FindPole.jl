


function GoStep(omgval::ComplexF64,
                IMat::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                FHT::FiniteHilbertTransform.FHTtype,
                tabaMcoef::Array{Float64,4},
                tabωminωmax::Matrix{Float64},
                Parameters::ResponseParameters)

    # fill the M matrix
    tabM!(omgval,MMat,tabaMcoef,tabωminωmax,FHT,Parameters)

    # compute the determinant of I-M
    detXival = detXi(IMat,MMat)

    # now compute a Jacobian via finite differences
    riomgval = omgval + 1.e-5
    upomgval = omgval + 1.e-5im

    tabM!(riomgval,MMat,tabaMcoef,tabωminωmax,FHT,Parameters)
    detXivalri = detXi(IMat,MMat)

    tabM!(upomgval,MMat,tabaMcoef,tabωminωmax,FHT,Parameters)
    detXivalup = detXi(IMat,MMat)

    dXirir = real(detXivalri-detXival)/1.e-5
    dXirii = imag(detXivalri-detXival)/1.e-5
    dXiupr = real(detXivalup-detXival)/1.e-5
    dXiupi = imag(detXivalup-detXival)/1.e-5

    # jacobian = [dXirir dXiupr ; dXirii dXiupi]
    # curpos = [real(detXival);imag(detXival)]
    # increment = jacobian \ (-curpos)

    increment1, increment2 = OrbitalElements.inverse2Dlinear(dXirir,dXiupr,dXirii,dXiupi,-real(detXival),-imag(detXival))

    # omgval += increment[1] + increment[2]*1im
    omgval += increment1 + increment2*1im
    if Parameters.VERBOSE > 0
        println("detXi=$detXival")
        println("omega=$omgval")
    end

    return omgval,detXival
end


function FindDeterminantZero(startingomg::ComplexF64,
                             IMat::Array{ComplexF64,2},MMat::Array{ComplexF64,2},
                             FHT::FiniteHilbertTransform.FHTtype,
                             tabaMcoef::Array{Float64,4},
                             tabωminωmax::Matrix{Float64},
                             Parameters::ResponseParameters;
                             TOL::Float64=1.e-12)

    # initial values
    detXival = 1.0
    omgval = startingomg

    while abs(detXival)^2 > TOL

        omgval,detXival = GoStep(omgval,IMat,MMat,FHT,tabaMcoef,tabωminωmax,Parameters)

        if abs(detXival) > 1.0
            break
        end
    end

    # what happens when this goes wrong?
    return omgval
end


function FindPole(startingomg::ComplexF64,
                  FHT::FiniteHilbertTransform.FHTtype,
                  Parameters::ResponseParameters,
                  TOL::Float64=1.e-12)

    # Preparinng computations of the response matrices
    MMat, tabaMcoef, tabωminωmax = PrepareM(Parameters)
    IMat = makeIMat(Parameters.nradial)

    bestomg = FindDeterminantZero(startingomg,IMat,MMat,FHT,tabaMcoef,tabωminωmax,Parameters,TOL=TOL)

    (Parameters.VERBOSE >= 0) && println("Best O for n1max=$(Parameters.n1max),nradial=$(Parameters.nradial) == $bestomg")

    return bestomg
end
