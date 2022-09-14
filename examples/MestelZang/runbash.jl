##################################################
#
#   Parsing of the command-line arguments
#
##################################################
using ArgParse # To parse command-line arguments
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    #####
    #   Basis
    #####
    "--G"
    help = "Gravitational constant"
    arg_type = Float64
    default = 1.0
    "--basisname"
    help = "Basis name"
    arg_type = String
    default = "CluttonBrock"
    "--rb"
    help = "Basis radius"
    arg_type = Float64
    default = 1.0
    "--kKA"
    help = "Basis parameter (only for Kalnajs basis)"
    arg_type = Int64
    default = 7
    "--lharmonic"
    help = "Considered harmonic number"
    arg_type = Int64
    default = 2
    "--nmax"
    help = "Number of radial basis elements"
    arg_type = Int64
    default = 100
    #####
    #   Potential
    #####
    "--modelname"
    help = "Model name"
    arg_type = String
    default = "Mestel"
    "--R0"
    help = "Characteristic length"
    arg_type = Float64
    default = 20.0
    "--V0"
    help = "Characteristic circular velocity"
    arg_type = Float64
    default = 1.0
    #####
    #   DF
    #####
    "--DFname"
    help = "DF name"
    arg_type = String
    default = "Zang"
    "--q0"
    help = "Characteristic temperature"
    arg_type = Float64
    default = 6.
    "--Rin"
    help = "Inner tapering radius"
    arg_type = Float64
    default = 1.0
    "--Rout"
    help = "Outer tapering radius"
    arg_type = Float64
    default = 11.5
    "--Rmax"
    help = "Maximal radius"
    arg_type = Float64
    default = 20.0
    "--xi"
    help = "Active fraction"
    arg_type = Float64
    default = 1.0
    "--mu"
    help = "Outer tapering exponant"
    arg_type = Int64
    default = 5
    "--nu"
    help = "Inner tapering exponant"
    arg_type = Int64
    default = 6
    #####
    #   Integration
    #####
    "--rmin"
    help = "Minimal considered radius"
    arg_type = Float64
    default = 0.1
    "--rmax"
    help = "Maximal considered radius"
    arg_type = Float64
    default = 20.0
    "--Ku"
    help = "Number of Legendre points"
    arg_type = Int64
    default = 100
    "--Kv"
    help = "Number of v points"
    arg_type = Int64
    default = 100
    "--Kw"
    help = "Number of w integration points"
    arg_type = Int64
    default = 100
    "--n1max"
    help = "Maximal radial resonance number"
    arg_type = Int64
    default = 1
    #####
    #   Complex plane image
    #####
    "--omgmin"
    help = "Minimal considered pattern speed"
    arg_type = Float64
    default = 0.5
    "--omgmax"
    help = "Maximal considered pattern speed"
    arg_type = Float64
    default = 1.5
    "--etamin"
    help = "Minimal considered growth rate"
    arg_type = Float64
    default = -0.1
    "--etamax"
    help = "Maximal considered growth rate"
    arg_type = Float64
    default = 0.5
    "--nomg"
    help = "Number of pattern speed points"
    arg_type = Int64
    default = 100
    "--neta"
    help = "Number of growth rate points"
    arg_type = Int64
    default = 100
    #####
    #   Other parameters
    #####
    "--verbose"
    help = "verbose"
    arg_type = Int64
    default = 1
    "--overwrite"
    help = "Overwritting existing computation files"
    arg_type = Bool
    default = false
end
parsed_args = parse_args(tabargs)


##################################################
#
#   Model construction
#
##################################################
import OrbitalElements
import AstroBasis
import FiniteHilbertTransform

##############################
# Basis
##############################
const G  = parsed_args["G"]
const basisname = parsed_args["basisname"]
const rb    = parsed_args["rb"]
const lmax  = parsed_args["lharmonic"]
const nmax  = parsed_args["nmax"]

if basisname == "CluttonBrock"
    const basis = AstroBasis.CB72Basis_create(lmax=lmax,nmax=nmax,G=G,rb=rb)
elseif basisname == "Kalnajs"
    const kKA = parsed_args["kKA"]
    const basis = basis = AstroBasis.K76Basis_create(lmax=lmax,nmax=nmax,G=G,rb=rb,kKA=kKA)
else
    error("Unknown basis name.")
end

##############################
# Model Potential
##############################
const modelname = parsed_args["modelname"]

const R0, V0 = parsed_args["R0"], parsed_args["V0"]
const ψ(r::Float64)::Float64   = OrbitalElements.ψMestel(r,R0,V0)
const dψ(r::Float64)::Float64  = OrbitalElements.dψMestel(r,R0,V0)
const d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψMestel(r,R0,V0)
const d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψMestel(r,R0,V0)
const d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψMestel(r,R0,V0)
const Ω₀ = OrbitalElements.Ω₀Mestel(R0,V0)

##############################
# Outputs directories
##############################
const wmatdir="wmat/"*basisname*"/"
const gfuncdir="gfunc/"*basisname*"/"
const modedir = "xifunc/"*basisname*"/"

##############################
# Model DF
##############################
const q0 = try convert(Int64,parsed_args["q0"]) catch; parsed_args["q0] end
const σ0 = OrbitalElements.sigmar_Mestel_DF(R0,V0,q0)
const C0 = OrbitalElements.normC_Mestel_DF(G,R0,V0,q0)

# Tapering radii
const Rin, Rout, Rmax  = parsed_args["Rin"], parsed_args["Rout"], parsed_args["Rmax"] 
# Self-gravity fraction
const xi = parsed_args["xi"]
# Tapering exponants                
const mu, nu = parsed_args["mu"], parsed_args["nu"]

const dfname = parsed_args["DFname"]*"_q_"*string(q0)*"_xi_"*string(xi)*"_mu_"*string(mu)*"_nu_"*string(nu)

const ndFdJ(n1::Int64,n2::Int64,
        E::Float64,L::Float64,
        ndotOmega::Float64)::Float64   = OrbitalElements.mestel_Zang_ndDFdJ(n1,n2,E,L,ndotOmega;
                                                                            R0=R0,Rin=Rin,Rmax=Rmax,
                                                                            V0=V0,
                                                                            xi=xi,C=C0,
                                                                            q=q0,sigma=σ0,
                                                                            nu=nu,mu=mu)

##############################
# Legendre / Integration parameters
##############################
# Radii for frequency truncations
const rmin, rmax = parsed_args["rmin"], parsed_args["rmax"]

const Ku, Kv, Kw = parsed_args["Ku"], parsed_args["Kv"], parsed_args["Kw"]
const FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

const lharmonic = lmax
const n1max = parsed_args["n1max"]

##################################################
#
#   Computations
#
##################################################
import CallAResponse

const verbose   = parsed_args["verbose"]
const overwrite = parsed_args["overwrite"]

# call the function to construct W matrices
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,wmatdir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,rmin,rmax,VERBOSE=verbose)

# call the function to compute G(u) functions
CallAResponse.RunGfunc(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,wmatdir,gfuncdir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,dfname,rmin,rmax,VERBOSE=verbose,OVERWRITE=overwrite)

# construct a grid of frequencies to probe
nbω0, nbη       = parsed_args["nomg"], parsed_args["neta"]
ω0min, ω0max    = parsed_args["omgmin"], parsed_args["omgmax"]
ηmin, ηmax      = parsed_args["etamin"], parsed_args["etamax"]

tabω = CallAResponse.gridomega(ω0min,ω0max,nbω0,ηmin,ηmax,nbη)
CallAResponse.RunM(tabω,dψ,d2ψ,gfuncdir,modedir,FHT,Kv,Kw,basis,lharmonic,n1max,Ω₀,modelname,dfname,rmin,rmax,VERBOSE=verbose,OVERWRITE=overwrite)
