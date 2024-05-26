"""

VERBOSE flag rules
0 : quiet
1 : progress tracking
2 : standard debug
4 : detailed timing outputs

"""

"""
    RunMatrices(ωlist::Array{ComplexF64}, FHT::FiniteHilbertTransform.AbstractFHT, params::LinearParameters)

Compute the response matrices for a given list of frequencies and parameters.

# Arguments
- `ωlist::Array{ComplexF64}`: An array of complex frequencies.
- `FHT::FiniteHilbertTransform.AbstractFHT`: An instance of a finite Hilbert transform.
- `params::LinearParameters`: A structure containing the linear parameters required for the computation.

# Returns
- `tabRMreal::Array{Float64, 3}`: A 3D array containing the real parts of the response matrices for each frequency.
- `tabRMimag::Array{Float64, 3}`: A 3D array containing the imaginary parts of the response matrices for each frequency.

# Description
This function prepares the computation of response matrices for a given list of complex frequencies (`ωlist`). It initializes necessary data structures, iterates over each frequency to compute the response matrix, and stores the real and imaginary parts of these matrices. The results are then saved to an HDF5 file.

# Details
- The function uses multithreading to parallelize the computation across available threads.
- The computed matrices are saved in an HDF5 file, with both real and imaginary parts stored separately.
- The function provides optional timing information for the second frequency computation if verbosity (`params.VERBOSE`) is enabled.

# Examples
```julia
# Example usage
ωlist = [1.0 + 0.1im, 2.0 + 0.2im, 3.0 + 0.3im]
FHT = FiniteHilbertTransform()
params = LinearParameters(nradial=3, VERBOSE=1)
tabRMreal, tabRMimag = RunMatrices(ωlist, FHT, params)
```
"""
function RunMatrices(ωlist::Array{ComplexF64},
    FHT::FiniteHilbertTransform.AbstractFHT,
    params::LinearParameters)

    # Validate inputs
    if isempty(ωlist)
        throw(ArgumentError("ωlist cannot be empty"))
    end
    if params.nradial <= 0
        throw(ArgumentError("params.nradial must be positive"))
    end

    # Preparing computations of the response matrices
    try
        tabMlist, tabaMcoef, tabωminωmax, FHTlist = PrepareM(Threads.nthreads(), FHT, params)
    catch e
        println("Error preparing matrices: ", e)
        return nothing
    end

    # Number of frequency values
    nω = length(ωlist)
    # Allocate containers for the matrices
    nradial = params.nradial
    tabRMreal = zeros(Float64, nradial, nradial, nω)
    tabRMimag = zeros(Float64, nradial, nradial, nω)

    (params.VERBOSE > 0) && println("LinearResponse.Xi.RunMatrices: computing $nω frequency values.")

    # Loop through all frequencies using multithreading with dynamic scheduling
    Threads.@threads for i = 1:nω
        k = Threads.threadid()

        try
            # Time the second frequency computation if verbosity is enabled
            if (i == 2) && (params.VERBOSE > 0)
                @time tabM!(ωlist[i], tabMlist[k], tabaMcoef, tabωminωmax, FHTlist[k], params)
            else
                tabM!(ωlist[i], tabMlist[k], tabaMcoef, tabωminωmax, FHTlist[k], params)
            end

            # Store real and imaginary parts of the matrix
            for q = 1:nradial
                for p = 1:nradial
                    tabRMreal[p, q, i] = real(tabMlist[k][p, q])
                    tabRMimag[p, q, i] = imag(tabMlist[k][p, q])
                end
            end
        catch e
            # what kinds of diagnostics would be helpful here for failure modes?
            (params.VERBOSE > 0) && println("Error in thread $k while processing frequency $i: ", e)
        end
    end

    # Save results to an HDF5 file
    try
        h5open(MatFilename(params), "w") do file
            # ω grid
            write(file, "omega", real(ωlist))
            write(file, "eta", imag(ωlist))
            # Matrices
            write(file, "MatricesR", tabRMreal)
            write(file, "MatricesI", tabRMimag)
            # Parameters
            WriteParameters(file, params)
        end
    catch e
        println("Error writing to HDF5 file: ", e)
    end

    return tabRMreal, tabRMimag
end
