


struct structPotentialtype

    modelname::String         # potential name (default "potential")
    G::Float64
    M::Float64
    bc::Float64
    Ω₀::Float64

    ψ::Function
    dψ::Function
    d2ψ::Function
    d3ψ::Function
    d4ψ::Function

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
