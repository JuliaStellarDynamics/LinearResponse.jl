"""
TODO:
1. what to do with default basis and FHT?

"""

struct LinearParameters

    # Orbital Elements parameters
    Orbitalparams::OrbitalElements.OrbitalParameters

    # Basis parameters
    # dimension and nradial are repeated (for allocation issues - intensively called)
    dimension::Int64
    nradial::Int64
    Basisparams::Dict{String,Union{String,Number}}

    # integration parameters
    Ku::Int64
    Kv::Int64
    Kw::Int64

    VMAPN::Int64
    ADAPTIVEKW::Bool

    KuTruncation::Int64

    # Files names, directories and handling
    modelname::String
    dfname::String

    wmatdir::String
    gfuncdir::String
    axidir::String
    modedir::String

    OVERWRITE::Bool

    # Resonances parameters
    lharmonic::Int64
    n1max::Int64

    nbResVec::Int64
    tabResVec::Matrix{Int64} # Matrix{Float64} == Array{Float64,2}

    # Other parameters
    VERBOSE::Int64
end

"""
    LinearParameters(basis;Orbitalparams,Ku,Kv,Kw,VMAPN,ADAPTIVEKW,KuTruncation,modelname,dfname,wmatdir,gfuncdir,modedir,OVERWRITE,lharmonic,n1max,VERBOSE)
"""
function LinearParameters(basis::AstroBasis.AbstractAstroBasis;
                          Orbitalparams::OrbitalElements.OrbitalParameters=OrbitalElements.OrbitalParameters(),
                          Ku::Int64=200,Kv::Int64=200,Kw::Int64=200,
                          VMAPN::Int64=1,ADAPTIVEKW::Bool=false,KuTruncation::Int64=10000,
                          modelname::String="model",dfname::String="df",
                          wmatdir::String="",gfuncdir::String="",axidir::String="",modedir::String="",OVERWRITE::Bool=false,
                          lharmonic::Int64=2,n1max::Int64=10,
                          VERBOSE::Int64=0)

    # Basis parameters
    Basisparams = AstroBasis.GetParameters(basis)

    # Basis parameters
    dimension = AstroBasis.dimension(basis)#Basisparams["dimension"]
    nradial = Basisparams["nradial"]

    # Resonance vectors
    nbResVec, tabResVec = MakeTabResVec(lharmonic,n1max,dimension)

    return LinearParameters(Orbitalparams,dimension,nradial,Basisparams,
                            Ku,Kv,Kw,VMAPN,ADAPTIVEKW,KuTruncation,
                            modelname,dfname,
                            wmatdir,gfuncdir,axidir,modedir,OVERWRITE,
                            lharmonic,n1max,nbResVec,tabResVec,
                            VERBOSE)
end
