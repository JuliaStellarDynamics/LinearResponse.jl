

"""
    getCMatrix(lmax::Int)

Computes the C matrix for the spherical harmonics Ylm(π/2, 0) up to a given maximum angular momentum `lmax`.

The C matrix contains the pre-computed values of Ylm(π/2, 0) for all l and m, where l ranges from 0 to `lmax` and m ranges from -l to l. These values are computed using the symmetry properties of the spherical harmonics and are stored in a matrix for efficient access.

# Arguments
- `lmax::Int`: The maximum angular momentum value up to which the C matrix is computed.

# Returns
The C matrix, a 2-dimensional array of size (lmax+1) x (2lmax+1), where each element represents the value of Ylm(π/2, 0) for a specific l and m.

# Notes
- The index for l starts from 0, so l=0 corresponds to index 1 in the C matrix.

# Computing Ylm(pi/2,0)
# Rather than load any special function libraries
# we pre-compute the values of Ylm(pi/2,0) necessary up to \ell=lmax
# We exploit the symmetry properties and the fact that we only need l-0,1,2,... and m=-l,-l+2,...,l
# See notes for details

BE VERY CAREFUL: l=0 corresponds to index=1

# Example

"""
function getCMatrix(lmax::Int64)
    CMatrix = zeros(Float64,lmax+1,lmax+1)
    # Making sure that the array is initially filled with zeros

    # !! ATTENTION, size of the array is (lmax+1) in both directions
    #for m=1:(lmax+1), l=1:(lmax+1)
    #    CMatrix[l,m] = 0.0
    #end

    # Ylm(pi/2,0) for (l,m)=(0,0)
    # the HARMONIC INDEX STARTS AT \ell=0
    CMatrix[1,1] = 1/sqrt(4pi)

    # First loop to fill in the diagonal elements
    for l=0:(lmax-1) #!! ATTENTION TO THE BOUNDS OF THE LOOP
        CMatrix[(l+1)+1,(l+1)+1] = - sqrt((2l+3.0)/(2l+2.0))*CMatrix[l+1,l+1]
    end

    # Second loop to fill in the remaining terms
    for m=0:lmax # We perform the loop at fixed m
        for l=m:2:(lmax-2) # At fixed m, we scan the \ell by step of 2 !! ATTENTION TO THE BOUNDS OF THE LOOP
            CMatrix[(l+2)+1,m+1] = - sqrt((2l+5.0)/(2l+1.0))*sqrt(((l-m+1.0)*(l+m+1.0))/((l-m+2.0)*(l+m+2.0)))*CMatrix[l+1,m+1]
        end
    end

    return CMatrix
end

"""
CYlm(CMatrix::Matrix{Float64}, l::Int64, m::Int64)

Compute the spherical harmonics prefactor for a given l and m.

This function returns the spherical harmonics prefactor based on the input CMatrix,
which is a matrix of precomputed values. The prefactor is calculated according to
the formula described in Fouvry et al. 2022, in Mean.jl.

Parameters
----------
CMatrix : Matrix{Float64}
    The matrix of precomputed values for the spherical harmonics.
l : Int64
    The degree of the spherical harmonics.
m : Int64
    The order of the spherical harmonics.

Returns
-------
Float64
    The spherical harmonics prefactor for the given l and m.


Wrapped function that returns the spherical harmonics prefactor
identical to Fouvry et al. 2022, in Mean.jl

Examples
--------


"""
function CYlm(CMatrix::Matrix{Float64},l::Int64,m::Int64)

    # the input l is shifted by +1 to retrieve the correct tabulated value
    if iseven(l)
        # Yl,m(pi/2,0) = Yl,-m(pi/2,0) for even l
        return CMatrix[l+1,abs(m)+1]
    else
        # Yl,m(pi/2) = - Yl,-m(pi/2,0) for odd l
        return CMatrix[l+1,abs(m)+1]*sign(m)
    end

end
