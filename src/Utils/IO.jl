

"""
load input configuration file
"""
function LoadConfiguration(inputfile::String)
    include(inputfile)
end

function CheckValidDirectory(dir::String)

    # check that this is in fact a directory
    if !isdir(dir)
        println("CallAResponse.IO.CheckValidDirectory: nonexistent directory $dir.")
        return -1
    end

    # check the trailing slash
    if dir[length(dir)]!='/'
        println("CallAResponse.IO.CheckValidDirectory: bad dir (needs trailing /) $dir")
        return -1
    end

    # finally, try opening a file to check write permissions
    try
        tst = open(dir*"tst.dat", "w")
        close(tst)
        rm(dir*"tst.dat")
    catch e
        println("CallAResponse.IO.CheckValidDirectory: cannot write test file to $dir")
        return -1
    end

    return 0

end

"""
check for existence of directories
"""
function CheckConfigurationDirectories(;wmatdir::String="",gfuncdir::String="",modedir::String="")

    # check for all specified directories (will succeed if no directory is specified, just might annoyingly print to wherever code is executed,)
    dircheckval = 0

    if !(wmatdir=="")
        dircheckval += CheckValidDirectory(wmatdir)
    end

    if !(gfuncdir=="")
        dircheckval += CheckValidDirectory(gfuncdir)
    end

    if !(modedir=="")
        dircheckval += CheckValidDirectory(modedir)
    end

    # could do something with the checks, but it's probably fine
    if dircheckval < 0
        return -1
    else
        return 0
    end

end


"""
    wmat_filename()

filename for a given wmat result
"""
function wmat_filename(wmatdir::String,
                        modelname::String,
                        lharmonic::Int64,n1::Int64,n2::Int64,
                        rb::Float64)

    return wmatdir*"wmat_"*modelname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rb)*".h5"
end

"""
    GFuncFilename()

filename for a given Gfunc result
"""
function GFuncFilename(gfuncdir::String,
                        modelname::String,
                        dfname::String,
                        lharmonic::Int64,n1::Int64,n2::Int64,
                        K_u::Int64,
                        rbasis::Float64)

    return gfuncdir*"Gfunc_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb"*string(rbasis)*".Ku"*string(K_u)*".h5"
end

"""
    mode_filename()

"""
function mode_filename(modedir::String,
                       modelname::String,
                       lharmonic::Int64,
                       n1max::Int64,K_u::Int64)

    return modedir*"ModeShape_"*modelname*"_l_"*string(lharmonic)*"_n1_"*string(n1max)*"."*string(K_u)*".h5"
end


"""
    det_filename()

"""
function DetFilename(modedir::String,
                      modelname::String,
                      dfname::String,
                      lharmonic::Int64,
                      n1max::Int64,K_u::Int64,rbasis::Float64)

    return modedir*"Determinant_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1max)*"_rb_"*string(rbasis)*"."*string(K_u)*".h5"
end

"""
wrtie the determinant array to a file
"""
function WriteDeterminant(detname::String,
                          tabomega::Array{Complex{Float64}},
                          tabdet::Array{Complex{Float64}})

    h5open(detname, "w") do file
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
                     K_u::Int64,
                     rbasis::Float64)

    return modedir*"TabAXi_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"_rb_"*string(rbasis)*"."*string(K_u)*".h5"
end


"""
    MFilename()

"""
function MFilename(modedir::String,
                   modelname::String,
                   dfname::String,
                   lharmonic::Int64,
                   K_u::Int64)

    return modedir*"TabM_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"."*string(K_u)*".h5"

end
