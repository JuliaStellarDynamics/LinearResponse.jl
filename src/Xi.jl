"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""


"""
    RunAXi(FHT,params)

function to make the decomposition coefficients "a" of the response matrix M

these values do not depend on the frequency being evaluated: which makes them good to table

is this struggling from having to pass around a gigantic array? what if we did more splitting?
"""
function RunAXi(FHT::FiniteHilbertTransform.AbstractFHT,
                params::LinearParameters)

    # check the directories + FHT values against the Parameters
    CheckDirectories(params.gfuncdir,params.axidir)
    CheckFHTCompatibility(FHT,params)
    
    # get relevant sizes
    nbResVec, tabResVec = params.nbResVec, params.tabResVec
    nradial  = params.nradial
    Ku       = FHT.Ku

    # allocate the workspace
    tabaMcoef = zeros(Float64,Ku,nradial,nradial)

    # allocate quadrature result and warnflag tables
    restab  = zeros(Float64,Ku)
    warntab = 0#zeros(Float64,Ku)


    # loop through all resonance vectors
    # make this parallel, but check that the loop is exited cleanly?
    # Threads.@threads for nres in 1:nbResVec
    for nres = 1:nbResVec

        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        (params.VERBOSE > 0) && println("LinearResponse.Xi.RunAXi: Starting on ($n1,$n2).")

        outputfilename = AxiFilename(n1,n2,params)
        # Check if the file already exist / has enough basis elements / overwritting imposed
        # false if no computation needed, then continue
        CheckFileNradial(outputfilename,params,"LinearResponse.Xi.RunAXi: ($n1,$n2) resonance") || continue

        # Read G(u) from the GFunc file for this resonance
        gfuncfilename  = GFuncFilename(n1,n2,params)
        tabGXi = h5read(gfuncfilename,"Gmat")

        for np = 1:nradial
            for nq = np:nradial
                # get the contribution
                restab,warn = FiniteHilbertTransform.GetaXi!(FHT,view(tabGXi,nq,np,:),restab,warntab)

                for k=1:Ku
                    # Warning if to many Inf or Nan values
                    #(warntab[k] > 3) && println("LinearResponse.Xi.RunAXi: NaN/Inf (warnflag=$(warntab[k])) values for (n1,n2)=($n1,$n2), (np,nq)=($np,$nq), and k=$k: $(restab[k]).")

                    # populate the symmetric matrix
                    tabaMcoef[k,nq,np] = restab[k] # Element (np,nq)
                    tabaMcoef[k,np,nq] = restab[k] # Element (nq,np). If np=nq overwrite (this is fine).
                end
            end
        end

        ωmin, ωmax = h5read(gfuncfilename,"omgmin"), h5read(gfuncfilename,"omgmax")
        # this is expensive enough to compute that we will want to save these
        # with the table fully constructed, loop back through to write after opening a file for the resonance
        h5open(outputfilename, "w") do outputfile
            # Mappings parameters
            write(outputfile, "omgmin",ωmin)
            write(outputfile, "omgmax",ωmax)
            # M_{n1,n2} decomposition coefficients
            write(outputfile,"aXi",tabaMcoef)
            # Parameters
            WriteParameters(outputfile,params)
        end
    end # resonance vector loop
end



"""
    StageAXi(params)
    
reads the decomposition's coefficients and extremal frequencies from HDF5 files
"""
function StageAXi(params::LinearParameters)

    # check the directories + basis and FHT values against the Parameters
    CheckDirectories(params.axidir)

    # get dimensions from the relevant tables
    nbResVec, tabResVec = params.nbResVec, params.tabResVec
    nradial  = params.nradial
    Ku       = params.Ku

    # allocate memory
    tabaMcoef = zeros(Float64,Ku,nradial,nradial,nbResVec)
    tabωminωmax = zeros(Float64,2,nbResVec)

    #Threads.@threads for nres=1:nbResVec # Loop over the resonances
    for nres = 1:nbResVec # Loop over the resonances
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres] # Current resonance (n1,n2)
        
        filename = AxiFilename(n1,n2,params)

        tmptabaMcoef = h5read(filename,"aXi")

        for np = 1:nradial
            for nq = 1:nradial
                for k=1:Ku
                    # fill in the (np,nq) value
                    tabaMcoef[k,nq,np,nres] = tmptabaMcoef[k,nq,np]
                end
            end
        end

        ωmin = h5read(filename,"omgmin")
        ωmax = h5read(filename,"omgmax")
        tabωminωmax[1,nres], tabωminωmax[2,nres] = ωmin, ωmax

    end # end resonance vector loop

    return tabaMcoef, tabωminωmax
end


"""

run the full linear response (Basis FT, G(u) computation and coefficient decomposition)
"""
function RunLinearResponse(model::OrbitalElements.Potential,
                            distributionfunction::DistributionFunction,
                            FHT::FiniteHilbertTransform.AbstractFHT,
                            basis::AstroBasis.AbstractAstroBasis,
                            params::LinearParameters)
    
    # call the function to construct W matrices
    RunWmat(model,FHT,basis,params)

    # call the function to compute G(u) functions
    RunGfunc(distributionfunction,FHT,params)

    # call the function to compute decomposition coefficients
    RunAXi(FHT,params)
end