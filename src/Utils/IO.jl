

"""
load input configuration file
"""
function LoadConfiguration(inputfile::String)
    include(inputfile)
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
function det_filename(modedir::String,
                      modelname::String,
                      dfname::String,
                      lharmonic::Int64,
                      n1max::Int64,K_u::Int64)

    return modedir*"Determinant_"*modelname*"_df_"*dfname*"_l_"*string(lharmonic)*"_n1_"*string(n1max)*"."*string(K_u)*".h5"
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
