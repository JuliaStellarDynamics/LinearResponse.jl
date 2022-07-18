
import AstroBasis
import OrbitalElements
using Printf


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



# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64       = OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr(r::Float64)::Float64    = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²(r::Float64)::Float64  = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Ω₀      =    OrbitalElements.isochrone_Omega0(bc,M,G)

a,e = 0.1, 0.1

# compute rperi and rapo
rp,ra = OrbitalElements.rpra_from_ae(a,e)
#println("rp=$rp ra=$ra")

# stuff in a hardwired version
rp=0.20317486071993035
ra=2.048615592190439

# put in some dummy values for testing: picking a resonance
n1 = -1
n2 = 2

Lval = OrbitalElements.L_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
Eval = OrbitalElements.E_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
println("E=$Eval L=$Lval")

Ω₁r,Ω₂r = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
ndotOmega = n1*Ω₁r + n2*Ω₂r

# compute dF/dJ: call out for value
valndFdJ  = ndFdJ(n1,n2,Eval,Lval,ndotOmega)
println("valndFdJ=$valndFdJ")
