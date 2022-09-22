


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
function WMatFilename(wmatdir::String,
                      modelname::String,
                      lharmonic::Int64,n1::Int64,n2::Int64,
                      rbasis::Float64,
                      Ku::Int64,Kv::Int64,Kw::Int64)

    return wmatdir*"wmat_"*modelname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*"_Kw_"*string(Kw)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(gfuncdir::String,
                        modelname::String,
                        dfname::String,
                        lharmonic::Int64,n1::Int64,n2::Int64,
                        rbasis::Float64,
                        Ku::Int64,Kv::Int64)

    return gfuncdir*"Gfunc_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*".h5"
end

"""
    mode_filename()

"""
function ModeFilename(modedir::String,
                       modelname::String,
                       dfname::String,
                       lharmonic::Int64,n1max::Int64,
                       rbasis::Float64,
                       Ku::Int64,Kv::Int64)

    return modedir*"ModeShape_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1max)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*".h5"
end


"""
    det_filename()

"""
function DetFilename(modedir::String,
                      modelname::String,
                      dfname::String,
                      lharmonic::Int64,n1max::Int64,
                      rbasis::Float64,
                      Ku::Int64,Kv::Int64)

    return modedir*"Determinant_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1max)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*".h5"
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
function AxiFilename(modedir::String,
                     modelname::String,
                     dfname::String,
                     lharmonic::Int64,
                     n1::Int64,n2::Int64,
                     rbasis::Float64,
                     Ku::Int64,Kv::Int64)

    return modedir*"TabAXi_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*".h5"
end


"""
    MFilename()

"""
function MFilename(modedir::String,
                   modelname::String,
                   dfname::String,
                   lharmonic::Int64,
                   rbasis::Float64,
                   Ku::Int64,Kv::Int64)

    return modedir*"TabM_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_rb_"*string(rbasis)*"_Ku_"*string(Ku)*"_Kv_"*string(Kv)*".h5"

end
