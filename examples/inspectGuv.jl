
using PyPlot
import OrbitalElements
import FiniteHilbertTransform
import AstroBasis
import CallAResponse

using HDF5



# set up the model
bc, M, G = 1.,1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)
rmin,rmax = 0.0,1.e5


# set up the sampling points
Ku = 200
FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)
tabu, Ku = FHT.tabu, FHT.Ku

# pick the resonance to consider
n1,n2 = -1,2

# set up model parameters
βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin=rmin,rmax=rmax)
αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)
ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

# set up the basis
rb        = 2.0   # the scale for the basis elements
lmax,nmax = 2,100  # usually lmax corresponds to the considered harmonics lharmonic
basis     = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
AstroBasis.fill_prefactors!(basis) # should this be in the creation by default?
ndim      = basis.dimension
nradial   = basis.nmax
basisFT   = CallAResponse.BasisFT_create(basis)


# set the wmatdir
wmatdir = "wmat/"
modelname = "PlummerE"
Kv = 200
Kw = 200
lharmonic=n2

# pick a u value to sample across v
tabu = [-0.65]
Wdata = CallAResponse.MakeWmatUV(ψ,dψ,d2ψ,d3ψ,βc,n1,n2,tabu,
                    Kv,lharmonic,basisFT,Ω₀=2.0,
                    rmin=rmin,
                    rmax=rmax,
                    Kw=Kw,
                    EDGE=0.01,
                    VERBOSE=1,
                    NINT=32)


# set some distribution function
function ndFdJ(n1::Int64,n2::Int64,
               E::Float64,L::Float64,
               ndotOmega::Float64;
               bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,Ra::Float64=1.0)

    return OrbitalElements.plummer_ROI_ndFdJ(n1,n2,E,L,ndotOmega,bc,M,astronomicalG,Ra)

end


# make the relevant prefactors
CMatrix = CallAResponse.getCMatrix(lharmonic)
pref    = -2.0*(2.0*pi)^(3)*CallAResponse.CYlm(CMatrix,lharmonic,n2)^(2)/(2.0*lharmonic+1.0)


# at fixed (np,nq), make the integral values G(u)
# although, you could modify this below to get more curves
np,nq = 1,1

uval = tabu[1] # only 1 value!
vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
deltav = (vmax-vmin) / (Kv)


for nq=1:1
    Jac1curve = zeros(Kv)
    Jac2curve = zeros(Kv)
    Jac3curve = zeros(Kv)
    fullcurve = zeros(Kv)
    basiscurve = zeros(Kv)
    for kvval=1:Kv
        vval = vmin + deltav*(kvval-0.5)

        Ω1,Ω2 = Wdata.tabΩ1Ω2[1,kvval,:]
        a,e = Wdata.tabAE[1,kvval,:]
        Eval,Lval = OrbitalElements.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=0.001)
        Jacαβ = OrbitalElements.JacαβToUV(n1,n2,ωmin,ωmax,vval)
        JacEL = tabJ[1,kvval]
        Jac1curve[kvval] = Jacαβ
        Jac2curve[kvval] = JacEL
        ndotΩ = n1*Ω1 + n2*Ω2
        Jac3curve[kvval] = ndFdJ(n1,n2,Eval,Lval,ndotΩ)
        basiscurve[kvval] = Wdata.tabW[np,1,kvval]*Wdata.tabW[nq,1,kvval]
        fullcurve[kvval] = Lval*(pref/Ω₀)*deltav*JacEL*Jacαβ*Wdata.tabW[np,1,kvval]*Wdata.tabW[nq,1,kvval]*ndFdJ(n1,n2,Eval,Lval,ndotΩ)*(1/Ω1)
    end

    # uncomment whichever line you want to see
    #plot(Wdata.tabUV[1,:,2],Jac1curve,color="black",linewidth=0.2)
    #plot(Wdata.tabUV[1,:,2],Jac2curve,color="black",linewidth=0.2)
    #plot(Wdata.tabUV[1,:,2],Jac3curve,color="black",linewidth=0.2)
    #plot(Wdata.tabUV[1,:,2],basiscurve,color="black",linewidth=0.2)
    #plot(Wdata.tabUV[1,:,2],fullcurve,color="black",linewidth=1.)
end

xlabel("v")
ylabel("quantity")
title("(n1,n2)=($n1,$n2), $(modelname)")
