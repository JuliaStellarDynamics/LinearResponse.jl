




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

    return tabomega

end
