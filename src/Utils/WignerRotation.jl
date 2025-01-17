using WignerD
using FastGaussQuadrature


# thesea are the functions for when g=sgn(x)
function Ifactor(lvalp::Int64,lvalq::Int64,n2::Int64,n3::Int64,N::Int64=10)
    CMatrix = getCMatrix(max(lvalp,lvalq))

    Yp = CYlm(CMatrix,lvalp,n2)
    Yq = CYlm(CMatrix,lvalq,n2)
    # Number of points for Gaussian quadrature
    # You can adjust N for higher accuracy
    # Get the Gaussian quadrature points and weights for the interval [0, pi]
    points, weights = gausslegendre(N)

    # Transform the points and weights to the interval [0, π]
    a, b = 0, π
    transformed_points = 0.5 * (b - a) * (points .+ 1) .+ a
    transformed_weights = 0.5 * (b - a) * weights

    # Define the integrand with the sin(beta) term
    function integrand(β,lvalp,lvalq,n2,n3)
        Rp = WignerD.wignerdjmn(lvalp, n2, n3, β)
        Rq = WignerD.wignerdjmn(lvalq, n2, n3, β)
        return sin(β) * Rp*Rq*sign(cos(β))
    end

    # Perform the quadrature sum
    result = sum(transformed_weights .* integrand.(transformed_points,lvalp,lvalq,n2,n3))

    return ((2π)^3) * Yp * Yq * result

end

function Jfactor(lvalp::Int64,lvalq::Int64,n2::Int64,n3::Int64,N::Int64=10)
    CMatrix = getCMatrix(max(lvalp,lvalq))

    Yp = CYlm(CMatrix,lvalp,n2)
    Yq = CYlm(CMatrix,lvalq,n2)

    # to check: is this the correct ordering of arguments for m,n?
    Rp = WignerD.wignerdjmn(lvalp, n2, n3, β)
    Rq = WignerD.wignerdjmn(lvalq, n2, n3, β)

    return 2n3 * ((2π)^3) * Yp * Yq * Rp * Rq

end


