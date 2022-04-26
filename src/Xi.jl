

function makeXiCoefficients!(tabaXi::Vector{Matrix{Vector{Float64}}},
                             tabResVec::Matrix{Int64},
                             tabnpnq::Matrix{Int64},
                             tabwGLquad::Vector{Float64},
                             tabPGLquad::Matrix{Float64},
                             tabINVcGLquad::Vector{Float64},
                             basedir::String="")

    # get relevant sizes
    K_u = size(tabaXi[1][1,1])[1]
    nb_npnq = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    for nresvec in 1:nbResVec
        n1,n2 = tabResVec[1,nresvec],tabResVec[2,nresvec]

        # open the resonance file
        filename = basedir*"gfunc/Gfunc_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5"
        file = h5open(filename,"r")

        for i_npnq=1:nb_npnq # Loop over the basis indices to consider
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)

            if (n1==0) & (n2==0)
                continue
            end
            if (n1==1) & (n2==-2)
                continue
            end

            # open the correct resonance vector
            # read in the correct G(u) function
            tabGXi = read(file,"GXinp"*string(np)*"nq"*string(nq))

            for k=0:(K_u-1) # Loop over the Legendre functions
                res = 0.0 # Initialisation of the result
                for i=1:K_u # Loop over the G-L nodes
                    w = tabwGLquad[i] # Current weight
                    G = tabGXi[i] # Current value of G[u_i]
                    P = tabPGLquad[k+1,i] # Current value of P_k. ATTENTION, to the shift of the array. ATTENTION, to the order of the arguments.
                    res += w*G*P # Update of the sum
                end
                res *= tabINVcGLquad[k+1] # Multiplying by the prefactor. ATTENTION, to the shift of the array
                #####
                # Adding the contribution to the result
                if (np == nq) # Case where we are on the diagonal
                    tabaXi[nresvec][np,nq][k+1] = res # Element (np,nq=np)
                else
                    tabaXi[nresvec][np,nq][k+1] = res # Element (np,nq)
                    tabaXi[nresvec][nq,np][k+1] = res # Element (nq,np)
                end
            end
        end
    end
end




"""
# Function that initialises a container tabXi
# to 0.0 + 0.0*im
"""
function tabXi_init!(tabXi::Array{Complex{Float64},2})
    for np=1:nradial # Loop over the index np
        for nq=1:nradial # Loop over the index nq
            tabXi[np,nq] = 0.0 + 0.0*im # Re-Initialising the array to 0.
        end
    end
end


"""
# Function that computes Xi[np,nq]
# for a given COMPLEX frequency omg
# in physical units, i.e. not rescaled by 1/Omega0.
# @IMPROVE -- The shape of the array could maybe be improved
"""
function tabXi!(omg::Complex{Float64},
                tabXi::Array{Complex{Float64},2},
                tabaXi::Vector{Matrix{Vector{Float64}}},
                tabResVec::Matrix{Int64},
                tabnpnq::Matrix{Int64},
                struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
                dpotential::Function,
                ddpotential::Function,
                LINEAR::String="Unstable",
                Omega0::Float64=1.0)

    # also needs
    # Omega0?

    nbResVec = size(tabResVec)[2]

    tabXi_init!(tabXi) # Initialising the array to 0.
    #####
    for nResVec=1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec] # Current resonance (n1,n2)

        #####
        omg_nodim = omg/Omega0 # Dimensionless frequency rescaled by Omega0
        varpi = OrbitalElements.get_varpi(omg_nodim,n1,n2,dpotential,ddpotential) # Getting the rescaled frequency. ATTENTION, the argument is the dimensionless frequency omg_nodim
        #####
        PerturbPlasma.get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR)
        tabDLeg = struct_tabLeg.tabDLeg # Name of the array where the D_k(w) are stored
        #####
        for i_npnq=1:nb_npnq # Loop over the basis indices to consider
            np, nq = tab_npnq[1,i_npnq], tab_npnq[2,i_npnq] # Current value of (np,nq)

            #####
            res = 0.0 + 0.0*im # Initialisation of the result
            #####
            for k=0:(K_u-1) # Loop over the Legendre functions
                res += tabaXi[nResVec][np,nq][k+1]*tabDLeg[k+1] # Adding a contribution. ATTENTION, to the shift of the array.
            end
            #####
            if (np == nq) # We are computing an element on the diagonal
                tabXi[np,nq] += res # Filling in the element (np,nq=np)
            else # We are computing a non-diagonal element
                tabXi[np,nq] += res # Filling in the element (np,nq)
                tabXi[nq,np] += res # Filling in the element (nq,np)
            end
        end
    end
end




# make an identity matrix
function makeIMat(nradial::Int64)
    IMat = zeros(Complex{Float64},nradial,nradial) # Static container for the identity matrix
    ##########
    for np=1:nradial # Loop over the radial elements to fill in the identity matrix
        IMat[np,np] = 1.0 + 0.0im # Creating the identity matrix
    end
    return IMat
end


function detXi(omg::Complex{Float64},
               tabXi::Array{Complex{Float64},2},
               tabaXi::Vector{Matrix{Vector{Float64}}},
               tabResVec::Matrix{Int64},
               tabnpnq::Matrix{Int64},
               struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
               dpotential::Function,
               ddpotential::Function,
               LINEAR::String="Unstable",
               Omega0::Float64=1.0)
    ####
    IMat = makeIMat(nradial)
    tabXi!(omg,tabXi,tabaXi,tabResVec,tab_npnq,struct_tabLeg,dpotential,ddpotential,LINEAR,Omega0)
    ####
    val = det(Symmetric(IMat-tabXi)) # Computing the determinant of (I-Xi). ATTENTION, we tell julia that the matrix is symmetric

    # only save the real portion
    return real(val) # Output
end
