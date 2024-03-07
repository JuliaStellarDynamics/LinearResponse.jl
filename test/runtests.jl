"""
for the radially-biased Plummer model: compute linear response theory

TEST against known values
"""

using AstroBasis
using FiniteHilbertTransform
using HDF5
using LinearResponse
using Test
using OrbitalElements


# Basis
G  = 1.
rb = 3.0
lmax,nradial = 2,5 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)

@testset "BasisTest" begin
    AstroBasis.tabUl!(basis,2,1.0)
    @test basis.tabUl[1] ≈ -0.3372710 atol=0.000001
    @test basis.tabUl[2] ≈ 0.5688238 atol=0.000001
end


# Model Potential
const modelname = "PlummerE"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.PlummerPotential()

@testset "PotentialTest" begin
    @test ψ(1.0,model) ≈ -0.707106 atol=0.000001
    @test dψ(1.0,model) ≈ 0.3535533 atol=0.000001
    @test d2ψ(1.0,model) ≈ -0.17677669 atol=0.000001
end

# Model Distribution Function
dfname = "roi0.75"

function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=0.75)

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# Linear Response integration parameters
Ku = 12    # number of Legendre integration sample points
Kv = 20    # number of allocations is directly proportional to this
Kw = 20    # number of allocations is insensitive to this (also time, largely)?


# Define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max = 1  # the Fiducial value is 10, but in the interest of a quick calculation, we limit ourselves to 1.

# output directories
wmatdir  = "./"
gfuncdir = "./"
modedir  = "./"

# Mode of response matrix computation
# Frequencies to probe
nOmega   = 280
Omegamin = -0.2
Omegamax = 0.2
nEta     = 280
Etamin   = -0.1
Etamax   = 0.1


VERBOSE   = 0
OVERWRITE = false
VMAPN     = 1
ADAPTIVEKW= false

OEparams = OrbitalElements.OrbitalParameters(EDGE=OrbitalElements.DEFAULT_EDGE,TOLECC=OrbitalElements.DEFAULT_TOLECC,TOLA=OrbitalElements.DEFAULT_TOLA,
                                             NINT=OrbitalElements.DEFAULT_NINT,
                                             da=OrbitalElements.DEFAULT_DA,de=OrbitalElements.DEFAULT_DE,
                                             ITERMAX=OrbitalElements.DEFAULT_ITERMAX,invε=OrbitalElements.DEFAULT_TOL)


Parameters = LinearResponse.LinearParameters(basis,Orbitalparams=OEparams,Ω₀=OrbitalElements.frequency_scale(model),Ku=Ku,Kv=Kv,Kw=Kw,
                                             modelname=modelname,dfname=dfname,
                                             wmatdir=wmatdir,gfuncdir=gfuncdir,modedir=modedir,axidir=modedir,
                                             lharmonic=lharmonic,n1max=n1max,
                                             VERBOSE=VERBOSE,OVERWRITE=OVERWRITE,
                                             VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW)



# package the Linear Response steps to compute M:
# 1. Call the function to construct W matrices
# 2. Run the G function calculation
# 3. Compute the matrix coefficients
println("Time to evaluate RunLinearResponse:")
@time LinearResponse.RunLinearResponse(model,ndFdJ,FHT,basis,Parameters)

# construct a grid of frequencies to probe
tabω = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

println("Time to constuct the matrics at all requested locations:")
@time tabRMreal, tabRMimag = LinearResponse.RunMatrices(tabω,FHT,Parameters)

println("Time to compute the determinant at all requested locations:")
@time tabdet = LinearResponse.RunDeterminant(tabω,FHT,Parameters)

# find a pole by using gradient descent
startingomg = 0.0 + 0.05im
@time bestomg,detval = LinearResponse.FindPole(startingomg,FHT,Parameters)
println("The zero-crossing frequency is $bestomg.")

@testset "ModeLocation" begin
    @test real(bestomg) ≈ 0.0 atol=1.e-10
    @test imag(bestomg) ≈ 0.043 atol=0.001
end

# for the minimum, go back and compute the mode shape
bestomgtest = 2.6047879816973595e-17 + 0.04334991803806192im
EV,EM = LinearResponse.ComputeModeTables(bestomgtest,FHT,Parameters)

modeRmin = 0.01
modeRmax = 15.0
nmode = 100
ModeRadius,ModePotentialShape,ModeDensityShape = LinearResponse.GetModeShape(basis,modeRmin,modeRmax,nmode,EM,Parameters)

@testset "ModeShape" begin
    @test ModeRadius[10] == 1.3727272727272728 # make sure we are testing the right location!
    @test real(ModePotentialShape[10]) ≈ -0.9597 atol=1.e-3
    @test imag(ModePotentialShape[10]) ≈ 0.0 atol=1.e-10
    @test real(ModeDensityShape[10]) ≈ 0.34074 atol=1.e-3
    @test imag(ModeDensityShape[10]) ≈ 0.0 atol=1.e-10
end