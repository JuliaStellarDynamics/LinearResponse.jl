
######################################
# Reusable files (Wmat, Gfunc and Axi)
######################################

"""
    WMatFilename()

filename for a given wmat result
"""
function WMatFilename(n1::Int64,n2::Int64,params::LinearParameters)

    bparams = params.Basisparams
    return params.wmatdir*"Wmat_"*params.modelname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(bparams["rb"])*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*"_Kw_"*string(params.Kw)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(n1::Int64,n2::Int64,params::LinearParameters)

    bparams = params.Basisparams
    return params.gfuncdir*"Gfunc_"*params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(bparams["rb"])*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

"""
    AxiFilename()

"""
function AxiFilename(n1::Int64,n2::Int64,params::LinearParameters)

    bparams = params.Basisparams
    return params.axidir*"Axi_"*params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(bparams["rb"])*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

##############################################
# Results files (Matrices, determinants, ....)
##############################################
"""
    OutputsEnd()
"""
function OutputsEnd(params::LinearParameters)

    bparams = params.Basisparams
    return params.modelname*"_df_"*params.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(params.n1max)*"_rb_"*string(bparams["rb"])*"_Ku_"*string(params.Ku)*"_Kv_"*string(params.Kv)*".h5"
end

"""
    ModeFilename()

"""
function ModeFilename(params::LinearParameters)

    return params.modedir*"ModeShape_"*OutputsEnd(params)
end


"""
    DetFilename()

"""
function DetFilename(params::LinearParameters)

    return params.modedir*"Determinant_"*OutputsEnd(params)
end

"""
    MFilename()

"""
function MatFilename(params::LinearParameters)

    return params.modedir*"Matrices_"*OutputsEnd(params)
end


#################
# Parameters dump
#################

"""
write all the parameters to a file
"""
function WriteParameters(filename::String,
                         params::LinearParameters,
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,params)
    end
end

function WriteParameters(file::HDF5.File,
                         params::LinearParameters)

    # Write orbital parameters
    OrbitalElements.WriteParameters(file,params.Orbitalparams)

    # Write basis parameters
    basisgroup = create_group(file,"BasisParameters")
    for (key,value) in params.Basisparams
        write(basisgroup,key,value)
    end

    # Write linear parameters
    group = create_group(file,"LinearParameters")
    for i = 1:fieldcount(LinearParameters)
        # If this field is a string, a number or a boolean, it is a parameter
        # (Prevent from dumping arrays)
        if fieldtype(LinearParameters,i) <: Union{String,Number,Bool}
            varname = string(fieldname(LinearParameters,i))
            write(group,varname,getfield(params,i))
        end
    end
end