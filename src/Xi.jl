"""

@IMPROVE: what is the best way to pass the format for reading in Gfunctions? (might want to consider a new Gfunc format)

"""


"""makeXiCoefficients!(tabaXi,tabResVec,tabnpnq,tabwGLquad,tabPGLquad,tabINVcGLquad,basedir)

function to make the coefficients in Xi
"""
function makeXiCoefficients!(tabaXi::Vector{Matrix{Vector{Float64}}},
                             tabResVec::Matrix{Int64},
                             tabnpnq::Matrix{Int64},
                             tabwGLquad::Vector{Float64},
                             tabPGLquad::Matrix{Float64},
                             tabINVcGLquad::Vector{Float64},
                             basedir::String="")

    # get relevant sizes
    K_u      = size(tabaXi[1][1,1])[1]
    nb_npnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]

    # loop through all resonance vectors
    for nresvec in 1:nbResVec
        n1,n2 = tabResVec[1,nresvec],tabResVec[2,nresvec]

        # open the resonance file
        filename = basedir*"gfunc/Gfunc_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5"
        file = h5open(filename,"r")

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tabnpnq[1,i_npnq], tabnpnq[2,i_npnq] # Current value of (np,nq)

            # open the correct resonance vector
            # read in the correct G(u) function
            tabGXi = read(file,"GXinp"*string(np)*"nq"*string(nq))

            # DEBUG nan values in tabGXi
            #sumG = sum(tabGXi)
            #println(println("n1=$n1, n2=$n2, np=$np, nq=$nq, sumG=$sumG."))

            # Loop over the Legendre functions
            for k=1:K_u

                res = 0.0 # Initialisation of the result

                for i=1:K_u # Loop over the G-L nodes
                    w = tabwGLquad[i] # Current weight
                    G = tabGXi[i] # Current value of G[u_i]
                    P = tabPGLquad[k,i] # Current value of P_k. ATTENTION, to the order of the arguments.
                    res += w*G*P # Update of the sum
                end

                res *= tabINVcGLquad[k] # Multiplying by the Legendre prefactor.

                # populate the symmetric matrix
                if (np == nq) # Case where we are on the diagonal
                    tabaXi[nresvec][np,nq][k] = res # Element (np,nq=np)
                else
                    tabaXi[nresvec][np,nq][k] = res # Element (np,nq)
                    tabaXi[nresvec][nq,np][k] = res # Element (nq,np)
                end

            end

            #h5write("xifunc/tabaXI_n1_"*string(n1)*"_n2_"*string(n2)*"_np_"*string(np)*"_nq_"*string(nq)*"."*string(K_u)*".h5","tabGXi",tabaXi[nresvec][np,nq])
        end
    end
end




"""tabXiInit!(tabXi)
Function that resets a container tabXi to 0.0 + 0.0*im
"""
function tabXiInit!(tabXi::Array{Complex{Float64},2})
    for np=1:nradial # Loop over the index np
        for nq=1:nradial # Loop over the index nq
            tabXi[np,nq] = 0.0 + 0.0*im # Re-Initialising the array to 0.
        end
    end
end


"""tabXi!(omg,tabXi,tabaXi,tabResVec,tabnpnq,struct_tabLeg,dpotential,ddpotential,LINEAR,Omega0)
Function that computes Xi[np,nq] for a given COMPLEX frequency omg in physical units, i.e. not (yet) rescaled by 1/Omega0.

@IMPROVE: The shape of the array could maybe be improved

See LinearTheory.jl for a similar version
"""
function tabXi!(omg::Complex{Float64},
                tabXi::Array{Complex{Float64},2},
                tabaXi::Vector{Matrix{Vector{Float64}}},
                tabResVec::Matrix{Int64},
                tabnpnq::Matrix{Int64},
                struct_tabLeg::PerturbPlasma.struct_tabLeg_type,
                dpotential::Function,
                ddpotential::Function,
                LINEAR::String="unstable",
                Omega0::Float64=1.0)

    # get dimensions from the relevant tables
    nb_npnq  = size(tabnpnq)[2]
    nbResVec = size(tabResVec)[2]
    K_u      = size(tabaXi[1][1,1])[1]

    tabXiInit!(tabXi) # Initialising the array to 0.
    #####
    for nResVec=1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nResVec], tabResVec[2,nResVec] # Current resonance (n1,n2)

        # Rescale to get the dimensionless frequency
        omg_nodim = omg/Omega0
        varpi = OrbitalElements.get_varpi(omg_nodim,n1,n2,dpotential,ddpotential) # Getting the rescaled frequency

        # get the Legendre integration values
        PerturbPlasma.get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR)
        tabDLeg = struct_tabLeg.tabDLeg # Name of the array where the D_k(w) are stored

        # Loop over the basis indices to consider
        for i_npnq=1:nb_npnq
            np, nq = tab_npnq[1,i_npnq], tab_npnq[2,i_npnq] # Current value of (np,nq)

            # Loop over the Legendre functions to add all contributions
            res = 0.0 + 0.0*im
            for k=1:K_u
                res += tabaXi[nResVec][np,nq][k]*tabDLeg[k]
            end

            # fill the full matrix
            if (np == nq) # We are computing an element on the diagonal
                tabXi[np,nq] += res # Filling in the element (np,nq=np)
            else # We are computing a non-diagonal element
                tabXi[np,nq] += res # Filling in the element (np,nq)
                tabXi[nq,np] += res # Filling in the element (nq,np)
            end

        end # basis index loop
        #h5write("xifunc/tabXI_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5","tabXi",tabXi)

    end # resonance loop

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
