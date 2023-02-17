

##############################
# Directories and files checks
##############################

function CheckDirectories(dirs...)

    # check that this is in fact a directory path ("/" at the end, "" return true)
    for dir in dirs
        if !isdirpath(dir)
            println("CallAResponse.IO.CheckDirectories: $dir is not an existing directory path.")
            return false
        end

        # finally, try opening a file to check write permissions
        try
            tst = open(dir*"tst.dat", "w")
            close(tst)
            rm(dir*"tst.dat")
        catch e
            println("CallAResponse.IO.CheckValidDirectory: cannot write test file to $dir")
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
                          Parameters::ResponseParameters,
                          preprint::String="")
        
    if isfile(filename)
        oldnradial = try h5read(filename,"ResponseParameters/nradial") catch; return true end
        if (Parameters.OVERWRITE == false) && (Parameters.nradial <= oldnradial)
            (Parameters.VERBOSE > 0) && println(preprint*" file already exists with higher nradial: no computation.")
            return false
        else
            (Parameters.VERBOSE > 0) && println(preprint*" file already exists (possibly with lower nradial) : recomputing and overwriting.")
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
function WMatFilename(n1::Int64,n2::Int64,Parameters::ResponseParameters)

    return Parameters.wmatdir*"Wmat_"*Parameters.modelname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*"_Kw_"*string(Parameters.Kw)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(n1::Int64,n2::Int64,Parameters::ResponseParameters)

    return Parameters.gfuncdir*"Gfunc_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end

"""
    AxiFilename()

"""
function AxiFilename(n1::Int64,n2::Int64,
                     Parameters::ResponseParameters)

    return Parameters.axidir*"Axi_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end

##############################################
# Results files (Matrices, determinants, ....)
##############################################
"""
    OutputsEnd()
"""
function OutputsEnd(Parameters::ResponseParameters)

    return Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(Parameters.n1max)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end

"""
    ModeFilename()

"""
function ModeFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"ModeShape_"*OutputsEnd(Parameters)
end


"""
    DetFilename()

"""
function DetFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"Determinant_"*OutputsEnd(Parameters)
end

"""
    MFilename()

"""
function MatFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"Matrices_"*OutputsEnd(Parameters)
end


#################
# Parameters dump
#################

"""
write all the parameters to a file
"""
function WriteParameters(filename::String,
                         Parameters::ResponseParameters,
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,Parameters)
    end
end

function WriteParameters(file::HDF5.File,
                         Parameters::ResponseParameters)

    group = create_group(file,"ResponseParameters")
    for i = 1:fieldcount(ResponseParameters)
        varname = string(fieldname(ResponseParameters,i))
        if (varname == "tabResVec")
            continue
        elseif (varname == "OEparams")
            OrbitalElements.WriteParameters(file,Parameters.OEparams)
        else
            try write(group,varname,getfield(Parameters,i)) catch; println("Unable to write parameter: "*varname) end
        end
    end
end