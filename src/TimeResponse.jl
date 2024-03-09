



"""
	inverseIMinusM!(MMat::Array{Array{Complex{Float64},2},1}, inverse::Array{Array{Complex{Float64},2},1})

Function filling the argument `inverse` with the inverse of the argument `matrix`. `matrix` is given as a list of square matrices, corresponding to the values of the response matrix at different times, and is considered as a single triangular block-Toeplitz matrix to compute the inverse. Since the inverse will also be a triangular block-Toeplitz matrix, the `inverse` also has the structure of a list of square matrices.

# Arguments
- `MMat`: list of square matrices, corresponding to the response matrix at different times.
- `inverse`: list of square matrices to fill.

# Output
None
"""
function inverseIMinusM!(MMat::Array{Array{Complex{Float64},2},1}, 
                         inverse::Array{Array{Complex{Float64},2},1},
                         DeltaT::Float64,
                         GridTimesSize::Int,
                         NBasisElements::Int)

    IMat = zeros(Complex{Float64},NBasisElements,NBasisElements)
    for np=1:NBasisElements
        IMat[np,np] = 1.0 + 0.0im 
    end
    inverse[1] = inv(IMat - DeltaT * MMat[1])
	for k=2:GridTimesSize
		matrixSum = zeros(Complex{Float64},NBasisElements,NBasisElements)
		for i=1:k-1
			matrixSum += - DeltaT * MMat[i + 1] * inverse[k - i]
		end
		inverse[k] = - inverse[1] * matrixSum
	end
end


"""
	response!(inverse::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})

Function filling the argument `response` with the product of the argument `inverse` minus the identity matrix with `perturber`. `inverse` is given as a list of square matrices, considered as a triangular block-Toeplitz matrix. 

# Arguments
- `inverse`: list of square matrices, corresponding to the matrix [(I-M)^-1] at different times.
- `perturber`: list of vectors corresponding to the perturber. The outer dimension concerns different time steps, while the inner dimension concerns the order of radial basis functions.
- `response`: list of vectors to fill.

# Output
None
"""
function response!(inverse::Array{Array{Complex{Float64},2},1}, 
                   perturber::Array{Array{Complex{Float64},1},1},
                   response::Array{Array{Complex{Float64},1},1},
                   GridTimesSize::Int,
                   NBasisElements::Int)

    IMat = zeros(Complex{Float64},NBasisElements,NBasisElements)
    for np=1:NBasisElements
        IMat[np,np] = 1.0 + 0.0im 
    end

	for k=1:GridTimesSize
		matrixSum = (inverse[1] - IMat) * perturber[k]
		for i=2:k
			matrixSum += inverse[i] * perturber[k - i + 1]
		end
		response[k] = matrixSum
	end
end

"""
	bareResponse!(matrix::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})

Function filling the argument `response` with the product of the argument `matrix` with `perturber`. `matrix` is given as a list of square matrices, considered as a triangular block-Toeplitz matrix. 

# Arguments
- `matrix`: list of square matrices, corresponding to the matrix [(I-M)^-1] at different times.
- `perturber`: list of vectors corresponding to the perturber. The outer dimension concerns different time steps, while the inner dimension concerns the order of radial basis functions.
- `response`: list of vectors to fill.

# Output
None
"""
function bareResponse!(matrix::Array{Array{Complex{Float64},2},1}, 
                       perturber::Array{Array{Complex{Float64},1},1},
                       response::Array{Array{Complex{Float64},1},1},
                       DeltaT::Float64,
                       GridTimesSize::Int)
	for k=1:GridTimesSize
		matrixSum = matrix[1] * perturber[k]
		for i=2:k
			matrixSum += DeltaT * matrix[i] * perturber[k - i + 1]
		end
		response[k] = matrixSum
	end
end



# Base cardinal sine and cosine are defined as
# Base.sinc(x) = sin(πx)/(πx)
# Base.cosc(x) = cos(πx)/x - sin(πx)/(πx^2)
# while we want
# sinc(t) = sin(t)/t
# cosc(t) = - sin(t)/t^2
mysinc(t::Float64) = Base.sinc(t/π)
mycosc(t::Float64) = Base.cosc(t/π)/π


function timeresponse!(res::Array{ComplexF64,1},t::Float64)

    # Initial values of PFT_k(t)
    val_0_PFT = 2im*mysinc(t) # Finite Fourier Transform of P_0
    val_1_PFT = - 2*mycosc(t)+0im # Finite Fourier Transform of P_1

    FiniteHilbertTransform.tabQLeg!(t+0im, val_0_PFT, val_1_PFT, res) # Computing the tabPFTLeg

    return res
end

function Getτ(
    tnodim::Number,
    ωmin::Float64,
    ωmax::Float64)

    return (ωmax - ωmin)*tnodim/2
end


function tabM!(t::Float64,
    tabM::AbstractMatrix{ComplexF64},
    tabaMcoef::Array{Float64,4},
    tabωminωmax::Matrix{Float64},
    fht::FiniteHilbertTransform.AbstractFHT,
    params::LinearParameters)

    # get dimensions from the relevant tables
    nbResVec, tabResVec = params.nbResVec, params.tabResVec
    KuTruncation = params.KuTruncation
    nradial  = params.nradial
    VERBOSE  = params.VERBOSE
    Ku       = fht.Ku

    # initialise the array to 0.
    fill!(tabM,0.0 + 0.0*im)

    # Rescale to get dimensionless frequency
    Ω₀ = params.Ω₀
    tnodim = Ω₀*t

    # loop over the resonances: no threading here because we parallelise over frequencies
    for nres = 1:nbResVec

        # get current resonance numbers (n1,n2)
        n1, n2 = tabResVec[1,nres], tabResVec[2,nres]

        # get ωmin and ωmax values
        ωmin, ωmax = tabωminωmax[1,nres], tabωminωmax[2,nres]

        # get the rescaled frequency
        τ = Getτ(tnodim,ωmin,ωmax)

        # get the integration values
        tabD = timeresponse!(fht.tabDLeg,τ)

        # Add the prefactors that are different from the frequency computations
        tabD .*= Ω₀ * (ωmax-ωmin)/2 * exp(-im*tnodim*(ωmax+ωmin)/2)

        # loop over the basis indices to consider
        for np = 1:nradial
            for nq = np:nradial

                res = 0.0 + 0.0*im

                # loop over the Legendre functions to add all contributions
                for k=1:Ku

                    # hard check for nans
                    @inbounds val = tabaMcoef[k,nq,np,nres]*tabD[k]
                    if !isnan(val)
                        res += val
                    else
                        (k==1) && (VERBOSE>1) && println("LinearResponse.Xi.tabM!: NaN found for n=($n1,$n2), npnq=($np,$nq), k=$k")
                    end


                end

                # fill the full M matrix:
                # as tab_npnq is the upper triangular matrix (with the diagonal),
                # we need to duplicate for symmetries

                # fill in the element (np,nq)
                @inbounds tabM[np,nq] += res

                # fill in the element (nq,np). @WARNING: added twice for diagonal elements np=nq.
                @inbounds tabM[nq,np] += res
            end
        end # basis index loop

        end # resonance loop

    # contributions were added twice for diagonal elements: reset
    for np=1:nradial
        @inbounds tabM[np,np] *= 0.5
    end
end