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

    nbResVec::Int64
    tabResVec::Matrix{Int64}

    tabnpnq::Matrix{Int64}

    KuTruncation::Int64

    VERBOSE::Int64
    OVERWRITE::Bool

    # Orbital Elements parameters
    Ω₀::Float64
    rmin::Float64
    rmax::Float64
    αmin::Float64
    αmax::Float64

    EDGE::Float64
    ELTOLECC::Float64

end


function ResponseParametersCreate(dψ::Function,d2ψ::Function,basis::AstroBasis.Basis_type;
                                  Ku::Int64=200,Kv::Int64=200,Kw::Int64=200,
                                  modelname::String="model",dfname::String="df",
                                  wmatdir::String="",gfuncdir::String="",modedir::String="",
                                  lharmonic::Int64=2,n1max::Int64=10,nradial::Int64=10,
                                  KuTruncation::Int64=10000,
                                  VERBOSE::Int64=0,OVERWRITE=false,
                                  Ω₀::Float64=1.0,rmin::Float64=1.e-5,rmax::Float64=1.e5,
                                  EDGE::Float64=0.01,ELTOLECC=0.001)

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,basis.ndim)

    # make the (np,nq) vectors that we need to evaluate
    tabnpnq = makeTabnpnq(nradial)

    # Frequency truncation parameters
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)

    return ResponseParameters(Ku,Kv,Kw,
                              modelname,dfname,
                              wmatdir,gfuncdir,modedir,
                              lharmonic,n1max,nradial,
                              nbResVec,tabResVec,
                              tabnpnq,
                              KuTruncation,
                              VERBOSE,OVERWRITE,
                              Ω₀,rmin,rmax,αmin,αmax,
                              EDGE,ELTOLECC)
end
