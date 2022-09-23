


struct structPotentialtype

    modelname::String         # potential name (default "potential")
    G::Float64                # gravitational constant (default 1.0)
    M::Float64                # model mass (default 1.0)
    bc::Float64               # model scale (default 1.0)
    Ω₀::Float64               # model frequency scale (default 1.0)

    ψ::Function               # model potential, to be queried at r::Float64 (default returns 0.0)
    dψ::Function              # model first derivative potential, to be queried at r::Float64 (default returns 0.0)
    d2ψ::Function             # model second derivative potential, to be queried at r::Float64 (default returns 0.0)
    d3ψ::Function             # model third derivative potential, to be queried at r::Float64 (default returns 0.0)
    d4ψ::Function             # model fourth derivative potential, to be queried at r::Float64 (default returns 0.0)

end


dummyfunc(r::Float64) = 0.

"""
    LegendreFHTcreate(Ku[name, dimension, lmax, nmax, G, rb])

Create a structLegendreFHTtype structure

"""
function structPotentialcreate(;modelname::String="potential",G::Float64=1.,bc::Float64=1.,M::Float64=1.,Ω₀::Float64=1.,
                                ψ::Function=dummyfunc,
                                dψ::Function=dummyfunc,
                                d2ψ::Function=dummyfunc,
                                d3ψ::Function=dummyfunc,
                                d4ψ::Function=dummyfunc)



    return structPotentialtype(modelname,G,M,bc,Ω₀,ψ,dψ,d2ψ,d3ψ,d4ψ)

end
