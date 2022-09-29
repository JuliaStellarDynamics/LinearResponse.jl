"""
TODO:
1. what to do with default basis and FHT?

"""

struct ResponseParameters

    # integration parameters
    Ku::Int64
    Kv::Int64
    Kw::Int64

    modelname::String
    dfname::String

    wmatdir::String
    gfuncdir::String
    modedir::String

    lharmonic::Int64
    n1max::Int64
    nradial::Int64
    ndim::Int64

    Ω₀::Float64
    rmin::Float64
    rmax::Float64

    basis::AstroBasis.Basis_type

    FHT::FiniteHilbertTransform.FHTtype

    nbResVec::Int64
    tabResVec::Matrix{Int64}

    tabnpnq::Matrix{Int64}

    KuTruncation::Int64

    VERBOSE::Int64
    EDGE::Float64
    ELTOLECC::Float64
    OVERWRITE::Boolean

end


function ResponseParametersCreate(basis,
                                  FHT;
                                  Ku=200,Kv=200,Kw=200,
                                  modelname="model",dfname="df",
                                  wmatdir="",gfuncdir="",modedir="",
                                  lharmonic=2,
                                  Ω₀=1.0,rmin=1.e-5,rmax=1.e5,
                                  n1max=10,
                                  nradial=10,
                                  ndim=3,
                                  KuTrunction=10000,
                                  VERBOSE=0,
                                  EDGE=0.01,ELTOLECC=0.001,OVERWRITE=false)

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    return ResponseParameters(Ku,Kv,Kw,
                              modelname,dfname,
                              wmatdir,gfuncdir,modedir,
                              lharmonic,n1max,nradial,ndim,
                              Ω₀,rmin,rmax,
                              basis,
                              FHT,
                              nbResVec,
                              tabResVec,
                              tabnpnq,
                              KuTrunction,
                              VERBOSE,
                              EDGE,ELTOLECC,OVERWRITE)
end
