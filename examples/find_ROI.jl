
import OrbitalElements
import PerturbPlasma
using HDF5
using LinearAlgebra

include("../src/Xi.jl")
include("../src/Basis.jl")         # define a list of the (np,nq) for which the response matrix must be computed
include("../src/Resonances.jl")    # for resonances helpers

basedir=""

#=
#const modelname = "isochrone"
const modelname = "isochroneE"

const bc, M, G = 1.,1.,1.
potential(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
=#

const modelname = "PlummerE"

const bc, M, G = 1.,1.,1.
potential(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
dpotential(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
ddpotential(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
Omega0 = OrbitalElements.plummer_Omega0(bc,M,G)

# define parameters: these should really be read in from somewhere else so there is no chance of confusion
K_u = 150
lmax=2
n1max=4
nradial = 10

# make the Legendre integration values
tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = PerturbPlasma.tabGLquad(K_u)

# make the resonance vectors
nbResVec = get_nbResVec(lmax,n1max)
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)

# make the (np,nq) vectors that we need to evaluate
tab_npnq = makeTabnpnq(nradial)

# initialise the array for memory saving
tabaXi = [[zeros(Float64,K_u) for np=1:nradial,nq=1:nradial] for ResVec=1:nbResVec]

# make the Xi coefficients
@time makeXiCoefficients!(tabaXi,tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,basedir)

#println(sum(tabaXi))

# once we have tabaXi, just go for it!
PARALLEL=true
struct_tabLeg = PerturbPlasma.initialize_struct_tabLeg(K_u,PARALLEL)


nbEta = 100 # Number of eta for which the ROI matrix is computed
Etamin = 0.001*Omega0 # Minimum eta value for which the ROI matrix is computed
Etamax = 0.030*Omega0 # Maximum eta value for which the ROI matrix is computed
deltaEta = (Etamax-Etamin)/(nbEta-1) # Getting the step distance of the array in eta
tabEta = collect(Etamin:deltaEta:Etamax) # Table of eta for which the response matrix is computed
tabdetXi = zeros(Float64,nbEta) # Table to store the value of det[I-Xi].


tabXi = zeros(Complex{Float64},nradial,nradial)

modeIndex = 1
#println("ModeIndex=$modeIndex")

# loop through values of eta
for iEta=1:nbEta
    omg = 0.0 + tabEta[iEta]*1.0im

    #if iEta==1
        tabdetXi[iEta] = detXi(omg,tabXi,tabaXi,tabResVec,tab_npnq,struct_tabLeg[1],dpotential,ddpotential,"unstable",Omega0)
    #end

    println(omg," ",tabdetXi[iEta])
    #if (abs(tabdetXi[iEta]) < abs(tabdetXi[modeIndex])) # We found an eta whose determinant is closer to 0
    #    modeIndex = iEta # Updating the index of the mode
    #end
end
#println("Eta unstable | ",tabEta[modeIndex])


##################################################
#print("Legendre projection | ") # Printing
#@time tabaXi!() # Computing the coefficients of the Legendre expansion

"""tabdetXi!()
Function that computes det[I-Xi] for all the considered frequencies

ATTENTION: we only keep the REAL part
"""
function tabdetXi!()
    Threads.@threads for iEta=1:nbEta # Loop over all the considered Eta
        thr = Threads.threadid() # ID of the current thread
        tabXi = struct_tabXi_parallel[thr].tabXi # Container for tabXi
        struct_tabLeg = struct_tabLeg_parallel[thr] # Struct container for the Legendre arrays
        #####
        eta = tabEta[iEta] # Current value of eta
        omg = 0.0 + im*eta # COMPLEX frequency considered
        #####
        val = detXi(omg,tabXi,struct_tabLeg) # Computing det[I-Xi] using the parallel containers
        #####
        tabdetXi[iEta] = real(val) # Filling in tabdetXi. ATTENTON, we only keep the real part.
    end
end
##################################################
print("Computation det[I-Xi] | ")
#@time tabdetXi!() # Computing all the det[I-Xi]
########################################
# Function that determines the eigenvector of the mode
########################################
tabEigenMode = zeros(Float64,nradial) # Constant container for the eigenmode
########################################
function tabEigenMode!()
    # First, we determine the index for which the determinant is the closest to zero
    modeIndex = 1 # Initial index
    for iEta=1:nbEta # Looping over all the Eta indices
        if (abs(tabdetXi[iEta]) < abs(tabdetXi[modeIndex])) # We found an eta whose determinant is closer to 0
            modeIndex = iEta # Updating the index of the mode
        end
    end
    println("Eta unstable | ",tabEta[modeIndex]) # Printing the eta of the unstable mode
    #####
    tabXi = struct_tabXi_serial.tabXi # Container for tabXi. ATTENTION, in serial.
    struct_tabLeg = struct_tabLeg_serial # Container for the Legendre arrays. ATTENTION, in serial.
    #####
    omgMode = 0.0 + tabEta[modeIndex]*im # COMPLEX frequency of the unstable mode
    #####
    tabXi!(omgMode,tabXi,struct_tabLeg) # Computing the response matrix of the unstable mode
    #####
    tabeigvals = eigvals(tabXi) # Table of the eigenvalues
    tabeigvecs = eigvecs(tabXi) # Table of the eigenvectors
    #####
    # Now we determine the eigenvalue that is the closest to 1.0
    indEig = 1 # Initial index
    for indrad=1:nradial # Loop over the number of basis elements
        if (abs(1.0-tabeigvals[indrad]) < abs(1.0-tabeigvals[indEig])) # We found an eigenvalue that is even closer to 1.0
            indEig = indrad # Updating the index of the eigenvalue
        end
    end
    ######
    # We fill in the eigenmode
    for np=1:nradial # Loop over the number of basis elements
        tabEigenMode[np] = real(tabeigvecs[np,indEig]) # Extracting the eigenvector !! ATTENTION, tabeigvecs is the matrix whose columns are eigenvectors !! ATTENTION, we put a real part
    end
end
#####
#tabEigenMode!() # Determining the eigenvector of the mode
