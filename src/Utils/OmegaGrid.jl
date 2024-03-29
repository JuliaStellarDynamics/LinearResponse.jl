


"""
make an identity matrix

@IMPROVE: where should this live?
"""
function makeIMat(nradial::Int64)
    IMat = zeros(ComplexF64,nradial,nradial) # Static container for the identity matrix
    ##########
    for np=1:nradial # Loop over the radial elements to fill in the identity matrix
        IMat[np,np] = 1.0 + 0.0im # Creating the identity matrix
    end
    return IMat
end

"""
make multiple identity matrices

@IMPROVE: where should this live?
"""
function makeIMat(nradial::Int64,nthreads::Int64)
    IMat = makeIMat(nradial)
    IMatlist = [deepcopy(IMat) for k=1:nthreads]
    return IMatlist
end

"""
function to make a table of omega values for frequency calculations
"""
function gridomega(Omegamin::Float64,Omegamax::Float64,nOmega::Int64,
                   Etamin::Float64,Etamax::Float64,nEta::Int64)

    # compute the total number of (complex) frequencies to probe
    nomega = nOmega*nEta

    # initialise the blank table
    tabomega = zeros(ComplexF64,nomega)

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


"""
function to make a table of omega values for frequency calculations
"""
function lineEta(Omegamin::Float64,Omegamax::Float64,nOmega::Int64,
                 Etaval::Float64;
                 log::Bool=false)

    # compute the total number of (complex) frequencies to probe
    nomega = nOmega

    # initialise the blank table
    tabomega = zeros(ComplexF64,nomega)

    # table of omega for which the response matrix is computed
    tabOmega = collect(range(Omegamin,Omegamax,length=nOmega))

    # initialise a counter
    icount = 1

    # loop over the real part of the frequency
    for iOmega=1:nOmega

        # fill the current value of the complex frequency
        if log
            tabomega[icount] = 10^(tabOmega[iOmega]) + im*Etaval
        else
            tabomega[icount] = tabOmega[iOmega] + im*Etaval
        end

        # update the counter
        icount += 1

    end

    return tabomega

end



"""
function to make a table of omega values for frequency calculations
"""
function lineOmega(Omegaval::Float64,
                   Etamin::Float64,Etamax::Float64,nEta::Int64;
                   log::Bool=false)

    # compute the total number of (complex) frequencies to probe
    nomega = nEta

    # initialise the blank table
    tabomega = zeros(ComplexF64,nomega)

    # table of eta for which the response matrix is computed
    tabEta = collect(range(Etamin,Etamax,length=nEta))

    # initialise a counter
    icount = 1

    # loop over the complex part of the frequency
    for iEta=1:nEta

        # fill the current value of the complex frequency
        if log
            tabomega[icount] = Omegaval + im*10^(tabEta[iEta])
        else
            tabomega[icount] = Omegaval + im*tabEta[iEta]
        end

        # update the counter
        icount += 1
    end

    return tabomega

end
