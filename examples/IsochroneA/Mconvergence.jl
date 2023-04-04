

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
using HDF5

function RunMIsochroneConvergence(inputfile::String,
                                  omglist::Array{Complex{Float64}};
                                  VERBOSE::Int64=0)

    # need some sort of 'if' for whether this already exists
    include(inputfile)

    #####
    # Check directories names
    #####
    if !(isdir(gfuncdir) && isdir(modedir))
        error("LinearResponse.Xi.RunMIsochrone: gfuncdir or modedir not found.")
    end

    for n1val=1:20
        # calculate the number of resonance vectors
        nbResVec = LinearResponse.get_nbResVec(lharmonic,n1val,ndim)

        # fill in the array of resonance vectors (n1,n2)
        tabResVec = LinearResponse.maketabResVec(lharmonic,n1val,ndim)

        # get all weights
        tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(K_u)

        # make the (np,nq) vectors that we need to evaluate
        tab_npnq = LinearResponse.makeTabnpnq(nradial)

        # make the decomposition coefficients a_k
        LinearResponse.MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,OVERWRITE=false,modedir=modedir)

        # allocate structs for D_k(omega) computation
        struct_tabLeglist = FiniteHilbertTransform.struct_tabLeg_create(K_u)

        # allocate memory for the response matrices M and identity matrices
        tabM = zeros(Complex{Float64},nradial,nradial)

        # make identity matrix and copies
        IMat = LinearResponse.makeIMat(nradial)

        # load aXi values
        tabaMcoef = LinearResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic)
        println("LinearResponse.Xi.RunMIsochrone: tabaMcoef loaded.")

        # loop through all frequencies
        #omgval = 0.0143 - 0.00142im # min for n1max=10; no .1. label
        omgval = 0.01 - 0.00073im # min for n1max=16; .2. label
        LINEAR = "damped"

        println("LinearResponse.Xi.RunMIsochrone: Starting frequency analysis, using $LINEAR integration.")
        #println("LinearResponse.Xi.RunMIsochrone: computing $nomglist frequency values.")

        @time LinearResponse.tabMIsochrone!(omgval,tabM,tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist,nradial,LINEAR,Omega0)

        detXival = LinearResponse.detXi(IMat,tabM)
        println("o=$omgval,n1max=$n1val,det=$detXival")

        h5open(modedir*"Mtab_n1max_"*string(n1val)*".2.hdf5", "w") do file
            write(file,"realM",real(tabM))
            write(file,"imagM",imag(tabM))
        end
        #println("n1max=$n1val,tabMlistval=$(tabMlist[1][nptest,nqtest])")

    end
        #WriteDeterminant(det_filename(modedir,modelname,dfname,lharmonic,n1max,K_u),omglist,tabdetXi)

        #return tabdetXi
end


inputfile = "ModelParamIsochrone_damped.jl"

include(inputfile)
tabomega = LinearResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

RunMIsochroneConvergence(inputfile,tabomega)
