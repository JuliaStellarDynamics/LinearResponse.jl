import OrbitalElements
import PerturbPlasma
using HDF5

include("/Users/mpetersen/CodeHold/JuliaCallAResponse/src/CMatrix.jl")
include("/Users/mpetersen/CodeHold/JuliaCallAResponse/src/GFunc.jl")

basedir="/Users/mpetersen/CodeHold/JuliaCallAResponse/examples/"


# bring in Legendre integration prefactors
K_u = 200
tabuGLquadtmp,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)
tabuGLquad = reshape(tabuGLquadtmp,K_u,1)

bc, M, G = 1.,1. ,1.
potential   = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential  = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Omega0      = OrbitalElements.isochrone_Omega0(bc,M,G)


lharmonic = 2
CMatrix = getCMatrix(lharmonic)
# example call: CYlm(CMatrix,1,0)

# decide on a resonances
for n1 in [-2,-1,0,1,2]
    for n2 in [-2,0,2]
        if (n1==0) & (n2==0)
            continue
        end

        # block this one for now...
        if (n1==1) & (n2==-2)
            continue
        end

        println(n1," ",n2)

        # load a value of tabWmat, plus (a,e) values
        filename = basedir*"wmat/wmat_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*".h5"
        file = h5open(filename,"r")
        Wtab = read(file,"wmat")
        atab = read(file,"amat")
        etab = read(file,"emat")
        nradial,K_u,K_v = size(Wtab)

        # need to loop through all combos of np and nq to make the full matrix.
        h5open(basedir*"gfunc/Gfunc_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5", "w") do file
            for np = 1:nradial
                for nq = 1:nradial
                    tabGXi = make_Gu(potential,dpotential,ddpotential,n1,n2,np,nq,Wtab,atab,etab,tabuGLquad,K_v,nradial,lharmonic,Omega0,bc,M,G)
                    write(file, "GXinp"*string(np)*"nq"*string(nq),tabGXi)
                end
            end
        end
    end
end
