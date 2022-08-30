

import OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import CallAResponse
using HDF5

function RunNIsochroneConvergence(inputfile::String,
                                  omglist::Array{Complex{Float64}};
                                  VERBOSE::Int64=0)

    # need some sort of 'if' for whether this already exists
    include(inputfile)

    #####
    # Check directories names
    #####
    if !(isdir(gfuncdir) && isdir(modedir))
        error("CallAResponse.Xi.RunMIsochrone: gfuncdir or modedir not found.")
    end

    # pick some fiducial n1val
    n1val = 8

    for nradial=1:100
        # calculate the number of resonance vectors
        nbResVec = CallAResponse.get_nbResVec(lharmonic,n1val,ndim)

        # fill in the array of resonance vectors (n1,n2)
        tabResVec = CallAResponse.maketabResVec(lharmonic,n1val,ndim)

        # get all weights
        tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = FiniteHilbertTransform.tabGLquad(K_u)

        # make the (np,nq) vectors that we need to evaluate
        tab_npnq = CallAResponse.makeTabnpnq(nradial)

        # make the decomposition coefficients a_k
        CallAResponse.MakeaMCoefficients(tabResVec,tab_npnq,tabwGLquad,tabPGLquad,tabINVcGLquad,gfuncdir,modedir,modelname,dfname,lharmonic,nradial,VERBOSE=VERBOSE,OVERWRITE=false,rb=rb)

        # allocate structs for D_k(omega) computation
        struct_tabLeglist = FiniteHilbertTransform.struct_tabLeg_create(K_u)

        # allocate memory for the response matrices M and identity matrices
        tabM = zeros(Complex{Float64},nradial,nradial)

        # make identity matrix and copies
        IMat = CallAResponse.makeIMat(nradial)

        # load aXi values
        tabaMcoef = CallAResponse.StageaMcoef(tabResVec,tab_npnq,K_u,nradial,modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)
        println("CallAResponse.Xi.RunMIsochrone: tabaMcoef loaded.")

        # chop the matrix down here and re-specify nradial

        # make the (np,nq) vectors that we need to evaluate
        tab_npnq = CallAResponse.makeTabnpnq(nradial)

        # loop through all frequencies
        omgval = 0.0143 - 0.00142im # min for n1max=10; no .1. label
        #omgval = 0.01 - 0.00073im # min for n1max=16; .2. label

        #println("CallAResponse.Xi.RunMIsochrone: computing $nomglist frequency values.")

        @time CallAResponse.tabMIsochrone!(omgval,tabM,tabaMcoef,tabResVec,tab_npnq,struct_tabLeglist,nradial,Omega0)

        detXival = CallAResponse.detXi(IMat,tabM)
        println("o=$omgval,nradial=$nradial,det=$detXival")

        h5open(modedir*"Mtab_n1max_"*string(n1val)*"_nradial_"*string(nradial)*".2.hdf5", "w") do file
            write(file,"realM",real(tabM))
            write(file,"imagM",imag(tabM))
        end
        #println("n1max=$n1val,tabMlistval=$(tabMlist[1][nptest,nqtest])")

    end

end


inputfile = "ModelParamIsochrone_damped.jl"

include(inputfile)
tabomega = CallAResponse.gridomega(Omegamin,Omegamax,nOmega,Etamin,Etamax,nEta)

RunNIsochroneConvergence(inputfile,tabomega)
