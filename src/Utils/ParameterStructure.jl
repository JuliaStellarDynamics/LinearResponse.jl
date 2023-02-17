"""
TODO:
1. what to do with default basis and FHT?

"""

struct ResponseParameters

    # basis needs to be separate (for multiple copies)
    #nmax::Int64

    # integration parameters
    Ku::Int64
    Kv::Int64
    Kw::Int64

    # Vector{Float64} == Array{Float64,1}
    #tabu::Vector{Float64}

    modelname::String
    dfname::String

    wmatdir::String
    gfuncdir::String
    axidir::String
    modedir::String

    lharmonic::Int64
    n1max::Int64
    nradial::Int64

    # Matrix{Float64} == Array{Float64,2}
    nbResVec::Int64
    tabResVec::Matrix{Int64}

    KuTruncation::Int64

    VERBOSE::Int64
    OVERWRITE::Bool

    # Orbital Elements parameters
    OEparams::OrbitalElements.OrbitsParameters

    # all Basis parameters shoujld be copied here
    nmax::Int64
    rbasis::Float64
    ndim::Int64

    VMAPN::Int64
    ADAPTIVEKW::Bool

end


function ResponseParametersCreate(OEparams::OrbitalElements.OrbitsParameters;
                                  Ku::Int64=200,Kv::Int64=200,Kw::Int64=200,
                                  modelname::String="model",dfname::String="df",
                                  wmatdir::String="",gfuncdir::String="",axidir::String="",modedir::String="",
                                  lharmonic::Int64=2,n1max::Int64=10,nradial::Int64=10,
                                  KuTruncation::Int64=10000,
                                  VERBOSE::Int64=0,OVERWRITE::Bool=false,ndim::Int64=3,
                                  nmax::Int64,rbasis::Float64,VMAPN::Int64=2,ADAPTIVEKW::Bool=false)

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,ndim)

    return ResponseParameters(Ku,Kv,Kw,
                              modelname,dfname,
                              wmatdir,gfuncdir,axidir,modedir,
                              lharmonic,n1max,nradial,
                              nbResVec,tabResVec,
                              KuTruncation,
                              VERBOSE,OVERWRITE,
                              OEparams,
                              nmax,rbasis,ndim,
                              VMAPN,ADAPTIVEKW)
end
