

##############################
# Directories and files checks
##############################

function CheckDirectories(dirs...)

    # check that this is in fact a directory path ("/" at the end, "" return true)
    for dir in dirs
        if !isdirpath(dir)
            println("LinearResponse.IO.CheckDirectories: $dir is not an existing directory path.")
            return false
        end

        # finally, try opening a file to check write permissions
        try
            tst = open(dir*"tst.dat", "w")
            close(tst)
            rm(dir*"tst.dat")
        catch e
            println("LinearResponse.IO.CheckValidDirectory: cannot write test file to $dir")
            return false
        end
    end
    return true
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

######################################
# Reusable files (Wmat, Gfunc and Axi)
######################################

"""
    WMatFilename()

filename for a given wmat result
"""
function WMatFilename(n1::Int64,n2::Int64,params::LinearParameters=LinearParameters())

    return params.wmatdir*"Wmat_"*params.modelname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(params.rbasis)*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*"_Kw_"*string(params.Kw)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(n1::Int64,n2::Int64,params::LinearParameters=LinearParameters())

    return params.gfuncdir*"Gfunc_"*params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(params.rbasis)*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

"""
    AxiFilename()

"""
function AxiFilename(n1::Int64,n2::Int64,params::LinearParameters=LinearParameters())

    return params.axidir*"Axi_"*params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(params.rbasis)*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

##############################################
# Results files (Matrices, determinants, ....)
##############################################
"""
    OutputsEnd()
"""
function OutputsEnd(params::LinearParameters=LinearParameters())

    return params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(params.n1max)*"_rb_"*string(params.rbasis)*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

"""
    ModeFilename()

"""
function ModeFilename(params::LinearParameters=LinearParameters())

    return params.modedir*"ModeShape_"*OutputsEnd(params)
end


"""
    DetFilename()

"""
function DetFilename(params::LinearParameters=LinearParameters())

    return params.modedir*"Determinant_"*OutputsEnd(params)
end

"""
    MFilename()

"""
function MatFilename(params::LinearParameters=LinearParameters())

    return params.modedir*"Matrices_"*OutputsEnd(params)
end


#################
# Parameters dump
#################

"""
write all the parameters to a file
"""
function WriteParameters(filename::String,
                         params::LinearParameters=LinearParameters(),
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,params)
    end
end

function WriteParameters(file::HDF5.File,
                         params::LinearParameters=LinearParameters())

    group = create_group(file,"LinearParameters")
    for i = 1:fieldcount(LinearParameters)
        varname = string(fieldname(LinearParameters,i))
        if (varname == "tabResVec")
            continue
        elseif (varname == "Orbitalparams")
            OrbitalElements.WriteParameters(file,params.Orbitalparams)
        else
            try write(group,varname,getfield(params,i)) catch; println("Unable to write parameter: "*varname) end
        end
    end
end