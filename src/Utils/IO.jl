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
    gfunc_filename()

filename for a given Gfunc result
"""
function gfunc_filename(gfuncdir::String,
                        modelname::String,
                        lharmonic::Int64,n1::Int64,n2::Int64,
                        K_u::Int64)

    return gfuncdir*"Gfunc_"*modelname*"_l_"*string(lharmonic)*"_n1_"*string(n1)*"_n2_"*string(n2)*"."*string(K_u)*".h5"
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
