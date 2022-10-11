

function SetupDeterminantZero(FHT::FiniteHilbertTransform.FHTtype,
                              ndim::Int64,lharmonic::Int64,
                              n1max::Int64,nradial::Int64,
                              modedir::String,modelname::String,
                              dfname::String,
                              rb::Float64,Kv::Int64;
                              VERBOSE::Int64=0)

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # read the decomposition coefficients a_k (must be created elsewhere)
    tabaMcoef = StageaMcoef(tabResVec,tabnpnq,FHT.Ku,Kv,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)

    # allocate arrays for constructing the determinant
    MMat = zeros(Complex{Float64},nradial,nradial)
    IMat = makeIMat(nradial)

    return IMat,MMat,tabaMcoef,tabResVec,tabnpnq

end


function GoStep(omgval::Complex{Float64},
                IMat::Array{Complex{Float64},2},MMat::Array{Complex{Float64},2},
                FHT::FiniteHilbertTransform.FHTtype,
                tabaMcoef::Array{Float64,4},
                tabResVec::Matrix{Int64},tabnpnq::Matrix{Int64},
                nradial::Int64,dψ::Function,
                d2ψ::Function,Ω₀::Float64,
                rmin::Float64,rmax::Float64;
                VERBOSE::Int64=0,
                KuTruncation::Int64=10000)

    nomg = 1

    tabdetXi = zeros(Float64,nomg) # real part of the determinant
    tabmevXi = zeros(Float64,nomg) # minimal eigenvalue at each frequency

    # fill the M matrix
    tabM!(omgval,MMat,tabaMcoef,tabResVec,tabnpnq,FHT,dψ,d2ψ,nradial,Ω₀,rmin,rmax,VERBOSE=VERBOSE,KuTruncation=KuTruncation)

    # compute the determinant of I-M
    detXival = detXi(IMat,MMat)

    # now compute a Jacobian via finite differences
    riomgval = omgval + 1.e-5
    upomgval = omgval + 1.e-5im

    tabM!(riomgval,MMat,tabaMcoef,tabResVec,tabnpnq,FHT,dψ,d2ψ,nradial,Ω₀,rmin,rmax,VERBOSE=VERBOSE,KuTruncation=KuTruncation)
    detXivalri = detXi(IMat,MMat)

    tabM!(upomgval,MMat,tabaMcoef,tabResVec,tabnpnq,FHT,dψ,d2ψ,nradial,Ω₀,rmin,rmax,VERBOSE=VERBOSE,KuTruncation=KuTruncation)
    detXivalup = detXi(IMat,MMat)

    dXirir = real(detXivalri-detXival)/1.e-5
    dXirii = imag(detXivalri-detXival)/1.e-5
    dXiupr = real(detXivalup-detXival)/1.e-5
    dXiupi = imag(detXivalup-detXival)/1.e-5

    jacobian = [dXirir dXiupr ; dXirii dXiupi]
    curpos = [real(detXival);imag(detXival)]
    increment = jacobian \ (-curpos)

    omgval += increment[1] + increment[2]*1im
    if VERBOSE > 0
        println("detXi=$detXival")
        println("omega=$omgval")
    end

    return omgval,detXival

end




function FindDeterminantZero(startingomg::Complex{Float64},
                             IMat::Array{Complex{Float64},2},MMat::Array{Complex{Float64},2},
                             FHT::FiniteHilbertTransform.FHTtype,
                             tabaMcoef::Array{Float64,4},
                             tabResVec::Matrix{Int64},tabnpnq::Matrix{Int64},
                             nradial::Int64,dψ::Function,
                             d2ψ::Function,Ω₀::Float64,
                             rmin::Float64,rmax::Float64;
                             TOL::Float64=1.0e-16,VERBOSE::Int64=0,
                             KuTruncation::Int64=10000)

    # initial values
    detXival = 1.0
    omgval = startingomg

    while abs(detXival)^2 > TOL

        omgval,detXival = GoStep(omgval,IMat,MMat,FHT,tabaMcoef,tabResVec,tabnpnq,nradial,dψ,d2ψ,Ω₀,rmin,rmax,VERBOSE=VERBOSE,KuTruncation=KuTruncation)

        if abs(detXival) > 1.0
            break
        end
    end

    # what happens when this goes wrong?
    return omgval

end



function FindPole(startingomg::Complex{Float64},
                  FHT::FiniteHilbertTransform.FHTtype,
                  dψ::Function,d2ψ::Function,
                  Parameters::ResponseParameters)

    IMat,MMat,tabaMcoef,tabResVec,tabnpnq = SetupDeterminantZero(FHT,Parameters.ndim,Parameters.lharmonic,
                                                                 Parameters.n1max,Parameters.nradial,Parameters.modedir,
                                                                 Parameters.modelname,Parameters.dfname,
                                                                 Parameters.rb,Parameters.Kv,VERBOSE=Parameters.VERBOSE)


    bestomg = FindDeterminantZero(startingomg,
                                  IMat,MMat,
                                  FHT,
                                  tabaMcoef,
                                  tabResVec,tabnpnq,
                                  Parameters.nradial,
                                  dψ,d2ψ,Parameters.Ω₀,Parameters.rmin,Parameters.rmax,VERBOSE=Parameters.VERBOSE,KuTruncation=Parameters.KuTruncation)

    (Parameters.VERBOSE >= 0) && println("Best O for n1max=$(Parameters.n1max),nradial=$(Parameters.nradial) == $bestomg")

    return bestomg

end
