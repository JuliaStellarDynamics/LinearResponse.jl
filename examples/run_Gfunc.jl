import OrbitalElements
import PerturbPlasma
using HDF5

include("../src/CMatrix.jl")
include("../src/GFunc.jl")
#include("../src/GFuncIsochrone.jl")

include("../src/Resonances.jl")    # for resonances helpers

basedir=""

# basis parameters to make sure we grab the right one
const rb = 10.


# bring in Legendre integration prefactors: pull from WMat
K_u = 150
tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

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

# specify which harmonic you wish to probe
lharmonic = 2
CMatrix = getCMatrix(lharmonic)
# example call: CYlm(CMatrix,1,0)



"""
wrap the ndFdJ function you want

a generic ndFdJ should take
n1,n2,E,L,ndotOmega

any other parameters should be wrapped into the function as external constants

"""
function ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)
    return ROIndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)
    #return ISOndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG)
end


"""
# FROM MEAN.JL
# Function that returns the value of n.dF/dJ
# where the DF follows the normalisation convention int dxdv F = Mtot
# Arguments are:
# + (n1,n2)
# + (E,L)
# + ndotOmega = n1*Omega1 + n2*Omega2
# ATTENTION, specific to an isotropic DF
"""
function ISOndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    dFdE = OrbitalElements.isochrone_isotropic_dDFdE(E,bc,M,astronomicalG) # Current value of dF/E. ATTENTION, the DF is assumed to be isotropic
    res = ndotOmega*dFdE # Current value of n.dF/dJ. ATTENTION, the DF is assumed to be isotropic
    #####
    return res
end

"""
# FROM ROI.JL
# Function that computes n.dF/dJ for the ROI DF
# ATTENTION, we put n.dF/dJ != 0 only for 0 < Q < 1
# @IMPROVE -- ideally we should have rather reduced
# the (alpha,beta) domain of integration
"""
function ROIndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.)

    Q = OrbitalElements.isochrone_Q_ROI(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    dFdQ = OrbitalElements.isochrone_Saha_dDFdQ(Q,Ra,bc,M,astronomicalG) # Value of dF/dQ
    dQdE, dQdL = OrbitalElements.isochrone_dQdE_ROI(E,L,Ra,bc,M,astronomicalG), OrbitalElements.isochrone_dQdL_ROI(E,L,Ra,bc,M,astronomicalG) # Values of dQ/dE, dQ/dL
    #####
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL) # Value of n.dF/dJ

    return res
end



# decide on a resonances
lmax  = 2  # maximum harmonic
n1max = 4  # maximum number of radial resonances to consider

nbResVec = get_nbResVec(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)

println(nbResVec)

#Threads.@threads for i = 1:nbResVec
for i = 1:nbResVec
    n1,n2 = tabResVec[1,i],tabResVec[2,i]
    println(n1," ",n2)

    # load a value of tabWmat, plus (a,e) values
    #filename = basedir*"wmat/wmat_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*".h5"
    filename = basedir*"wmat/wmat_"*string(modelname)*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rb)*".h5"
    file = h5open(filename,"r")
    Wtab = read(file,"wmat")
    atab = read(file,"amat")
    etab = read(file,"emat")
    nradial,K_u,K_v = size(Wtab)
    println("nradial=$nradial,K_u=$K_u,K_v=$K_v")

    # Overall prefactor of the integrand. ATTENTION, to the minus sign.
    #pref = -2.0*(2.0*pi)^(3)*CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*lharmonic+1.0)

    # need to loop through all combos of np and nq to make the full matrix.
    h5open(basedir*"gfunc/Gfunc_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5", "w") do file
        for np = 1:nradial
            for nq = 1:nradial
                #@time tabGXi = makeGu(potential,dpotential,ddpotential,ndFdJ,n1,n2,np,nq,Wtab,atab,etab,tabuGLquad,K_v,nradial,lharmonic,pref,Omega0=Omega0,bc=bc,M=M,G=G)
                tabGXi = makeGu(potential,dpotential,ddpotential,ndFdJ,n1,n2,np,nq,Wtab,atab,etab,tabuGLquad,K_v,nradial,lharmonic,ndim=3,Omega0=Omega0,bc=bc,M=M,G=G)
                sumG = sum(tabGXi)
                if (np>-100) & (nq>-100)
                    if isnan(sumG)
                        println("NaN for n1=$n1, n2=$n2.")
                    else
                        #println("np=$np, nq=$nq, sumG=$sumG.")
                    end
                end
                write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)
            end
        end
    end
end
