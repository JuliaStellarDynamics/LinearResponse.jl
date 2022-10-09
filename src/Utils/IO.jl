


function CheckValidDirectory(dir::String)

    # check that this is in fact a directory
    if !isdir(dir)
        println("CallAResponse.IO.CheckValidDirectory: nonexistent directory $dir.")
        return false
    end

    # check the trailing slash
    if last(dir)!='/'
        println("CallAResponse.IO.CheckValidDirectory: bad dir (needs trailing /) $dir")
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

    return true
end

"""
check for existence of directories
"""
function CheckConfigurationDirectories(dirs::Array{String})

    # check for all specified directories (will succeed if no directory is specified, just might annoyingly print to wherever code is executed,)
    for dir in dirs
        ((dir=="") || CheckValidDirectory(dir)) || (return false)
    end
    return true
end


"""
    WMatFilename()

filename for a given wmat result
"""
function WMatFilename(n1::Int64,n2::Int64,Parameters::ResponseParameters)

    return Parameters.wmatdir*"wmat_"*Parameters.modelname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*"_Kw_"*string(Parameters.Kw)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(n1::Int64,n2::Int64,Parameters::ResponseParameters)

    return Parameters.gfuncdir*"Gfunc_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end

"""
    mode_filename()

"""
function ModeFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"ModeShape_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(Parameters.n1max)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end


"""
    det_filename()

"""
function DetFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"Determinant_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(Parameters.n1max)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end

"""
wrtie the determinant array to a file
"""
function WriteDeterminant(detname::String,
                          tabomega::Array{Complex{Float64}},
                          tabdet::Array{Complex{Float64}})

    h5open(detname, "w") do file
        #write(file,"nx",)
        write(file,"omega",real(tabomega))
        write(file,"eta",imag(tabomega))
        write(file,"det",tabdet)
    end
end

"""
    AxiFilename()

"""
function AxiFilename(n1::Int64,n2::Int64,
                     Parameters::ResponseParameters)

    return Parameters.modedir*"TabAXi_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"
end


"""
    MFilename()

"""
function MFilename(Parameters::ResponseParameters)

    return Parameters.modedir*"TabM_"*Parameters.modelname*"_df_"*Parameters.dfname*"_l_"*string(Parameters.lharmonic)*"_rb_"*string(Parameters.rbasis)*"_Ku_"*string(Parameters.Ku)*"_Kv_"*string(Parameters.Kv)*".h5"

end
