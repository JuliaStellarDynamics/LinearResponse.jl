


nOmega = 52#501 # Number of omega0 for which the response matrix is computed. ATTENTION, we choose this number to have a nice step distance
Omegamin = 0.0#*Omega0 # Minimum omega value for which the response matrix is computed. ATTENTION, we use the symmetry of the response matrix, so that we only search for positive real(omega).
Omegamax = 0.05#1*Omega0 # Maximum omega value for which the response matrix is computed

# Evaluate the dispersion function even in eta=0.0 using the `damped' evaluation,
# but this should not really change the plot that we are making in the lower-half of the complex plane.
nEta = 51#501 # Number of eta for which the response matrix is computed. ATTENTION, we choose this number to have a nice step distance
Etamin = -0.005#-0.06*Omega0#-0.005*Omega0 # Minimum eta value for which the response matrix is computed. ATTENTION, it is a negative number
Etamax = 0.0#-0.0*Omega0 # Maximum eta value for which the response matrix is computed. ATTENTION, it is a negative number
#####

"""
function to make a table of omega values for frequency calculations
"""
function gridomega(Omegamin::Float64,Omegamax::Float64,nOmega::Int64,
                   Etamin::Float64,Etamax::Float64,nEta::Int64)

    # compute the total number of (complex) frequencies to probe
    nomega = nOmega*nEta

    # initialise the blank table
    tabomega = zeros(Complex{Float64},nomega)

    # table of omega for which the response matrix is computed
    tabOmega = collect(range(Omegamin,Omegamax,length=nOmega))

    # table of eta for which the response matrix is computed
    tabEta = collect(range(Etamin,Etamax,length=nEta))


    # initialise a counter
    icount = 1

    # loop over the real part of the frequency
    for iOmega=1:nOmega

        # loop over the complex part of the frequency
        for iEta=1:nEta

            # fill the current value of the complex frequency
            tabomega[icount] = tabOmega[iOmega] + im*tabEta[iEta]

            # update the counter
            icount += 1
        end
    end

end
