

##############################
# Directories and files checks
##############################

function CheckDirectories(dirs...)

    # check that this is in fact a directory path ("/" at the end, "" return true)
    for dir in dirs
        if !isdirpath(dir)
            error("LinearResponse.IO.CheckDirectories: $dir is not an existing directory path.")
        end

        # finally, try opening a file to check write permissions
        try
            tst = open(dir*"tst.dat", "w")
            close(tst)
            rm(dir*"tst.dat")
        catch e
            error("LinearResponse.IO.CheckValidDirectory: cannot write test file to $dir")
        end
    end
end

"""
    Check if a file exist, has enough basis elements / if overwritting is expected
    Return true if computation needed, and false if not, i.e. if all following conditions are satisfied:
        - the file already exists,
        - overwritting is not mandatory (OVERWRITE == false)
        - the old file has been created with enough basis elements (nradial sufficiently high)
"""
function CheckFileNradial(filename::String,
                          params::LinearParameters=LinearParameters(),
                          preprint::String="")
        
    if isfile(filename)
        oldnradial = try h5read(filename,"LinearParameters/nradial") catch; return true end
        if (params.OVERWRITE == false) && (params.nradial <= oldnradial)
            (params.VERBOSE > 0) && println(preprint*" file already exists with higher nradial: no computation.")
            return false
        else
            (params.VERBOSE > 0) && println(preprint*" file already exists (possibly with lower nradial) : recomputing and overwriting.")
            return true
        end
    else 
        return true
    end
end

##############################
# FHT vs Parameters
##############################
function CheckFHTCompatibility(FHT::FiniteHilbertTransform.AbstractFHT,
                          params::LinearParameters)
        
    @assert (FHT.Ku == params.Ku) "Incompatible FHT and parameters."
end

##############################
# Basis vs Parameters
##############################
function CheckBasisCompatibility(basis::AstroBasis.AbstractAstroBasis,
                                 params::LinearParameters)

    compat = true
    for (key,value) in params.Basisparams
        if value != getfield(basis,Symbol(key))
            compat = false
            break
        end
    end
    @assert compat "Incompatible basis and parameters."
end

####################################
# Matrices shapes control
####################################
function CheckMShape(tabM::AbstractMatrix;
                             nradial::Int64)

    @assert size(tabM) == (nradial,nradial) "tabM shape not suitable."
end

function CheckaMcoefShape(tabaMcoef::Array{Float64,4};
                             nradial::Int64,Ku::Int64,nbResVec::Int64)

    k1, k2, k3, k4 = size(tabaMcoef)
    @assert ((k1 == Ku) && (k2 >= nradial) && (k3 >= nradial) && (k4 == nbResVec)) "tabaMcoef shape not suitable."
end