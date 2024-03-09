"""
for the radially-biased Isochrone model: compute linear response theory
"""

using OrbitalElements
using AstroBasis
using FiniteHilbertTransform
using LinearResponse
using HDF5

using Plots


# Basis
G  = 1.
rb = 2.0
lmax,nradial = 2,5 # Usually lmax corresponds to the considered harmonics lharmonic
basis = AstroBasis.CB73Basis(lmax=lmax, nradial=nradial,G=G,rb=rb)



# Model Potential
const modelname = "IsochroneE"
const bc, M = 1.,1. # G is defined above: must agree with basis!
model = OrbitalElements.NumericalIsochrone()

rmin = 0.0
rmax = Inf


dfname = "roi1.0"

function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.0)

    Q = OrbitalElements.isochroneQROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochroneSahadDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochronedQdEROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochronedQdLROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res

end


# Linear Response integration parameters
Ku = 50    # number of Legendre integration sample points
Kv = 30    # number of allocations is directly proportional to this
Kw = 20    # number of allocations is insensitive to this (also time, largely)?


# Define the helper for the Finite Hilbert Transform
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

lharmonic = lmax
n1max = 1 

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


VERBOSE   = 2
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
#@time LinearResponse.RunLinearResponse(model,ndFdJ,FHT,basis,Parameters)

@time LinearResponse.RunWmat(model,FHT,basis,Parameters)

# call the function to compute G(u) functions
@time LinearResponse.RunGfunc(ndFdJ,FHT,Parameters)

# call the function to compute decomposition coefficients
@time LinearResponse.RunAXi(FHT,Parameters)

# reload as an array
MMat, tabaMcoef, tabωminωmax = LinearResponse.PrepareM(Parameters)

# call the M matrix for each of the times and the basis elements
NBasisElements = nradial
GridTimesSize = 600
DeltaT = 0.5


# calculate M matrix with times
MMat = [zeros(Complex{Float64}, NBasisElements,NBasisElements) for _ in 1:GridTimesSize]

# the first matrix needs to be zero (as per Simon)
# loop through times
for t=2:GridTimesSize
    #MMat[t] = ones(Complex{Float64}, NBasisElements,NBasisElements)
    LinearResponse.tabM!((t-1)*DeltaT,MMat[t],tabaMcoef,tabωminωmax,FHT,Parameters)
    MMat[t] = real.(MMat[t])
end

println("Check reality")
println(MMat[290][1,1])
println(MMat[290][1,2])
println(MMat[290][1,3])
println(MMat[290][1,5])

# define the storage for the inverse
inverse = [zeros(Complex{Float64}, NBasisElements,NBasisElements) for _ in 1:GridTimesSize]

# invert the list of M matrices
LinearResponse.inverseIMinusM!(MMat,inverse,DeltaT,GridTimesSize,NBasisElements)

# define the basic perturber
onesPerturber = [exp(-(DeltaT * i)^2/40) .* ones(Complex{Float64}, NBasisElements) for i=1:GridTimesSize]

# compute the response given the perturber
response = [zeros(Complex{Float64}, NBasisElements) for _ in 1:GridTimesSize]
LinearResponse.response!(inverse,onesPerturber,response,GridTimesSize,NBasisElements)



# now, measure the growth rate of a single element of the response matrix
# response[time][basiselement]
BasisElement = 1 # select the basis element


timevals = 1:GridTimesSize
colors = cgrad(:GnBu_5,rev=true)

for BasisElement=1:NBasisElements
    timerun = zeros(GridTimesSize)
    for t=1:GridTimesSize
        timerun[t] = abs(real(response[t][BasisElement])) # @ATTENTION, is this real or imaginary? Or neither? absolute value?
    end
    if BasisElement==1
        plot(timevals,log.(timerun),linecolor=colors[BasisElement],label="p=$BasisElement")
        #plot(timevals,(timerun),linecolor=colors[BasisElement],label="p=$BasisElement")
    else        
        plot!(timevals,log.(timerun),linecolor=colors[BasisElement],label="p=$BasisElement")
        #plot!(timevals,(timerun),linecolor=colors[BasisElement],label="p=$BasisElement")
    end        
    # let's do finite difference to be simple, of some fairly late time
    DerivTime = 590
    finitediffslope = (log(timerun[DerivTime+1])-log(timerun[DerivTime-1]))/(2*DeltaT)
    println("Finite difference derivative for basis element $BasisElement at t=$DerivTime is γ=$finitediffslope")

end

# Add labels and title
xlabel!("time")
ylabel!("log amplitude")
savefig("ROItimedependent.png")


