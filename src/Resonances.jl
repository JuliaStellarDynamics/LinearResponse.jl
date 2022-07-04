"""
helper functions for making resonance containers

example:
lmax=2
n1max=2
nbResVec = get_nbResVec(lmax,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lmax
tabResVec = maketabResVec(lmax,n1max) # Filling in the array of resonance vectors (n1,n2)
println(tabResVec)

"""



"""
# Function that returns the total number
# of resonances bn=(n1,n2) to consider
# for the harmonics lmax
# There a few constraints to satisfy:
# + |n2| <= lmax
# + (lmax-n2) even
# + |n1| <= n1max
# + (n1,n2) = (0,0) does not contribute
# ATTENTION, the (n1,n2) are determined for l=lmax
"""
function get_nbResVec(lmax::Int64,n1max::Int64)
    count = 0 # Initialisation of the counter
    #####
    for n2=lmax:-2:0 # Loop over the index n2. ATTENTION, we go by step of 2, and only consider positive values
        #####
        n1bound = -n1max # Value of the bound for n1 if n2 > 0
        #####
        if (n2 == 0)
            n1bound = 1 # Value of the bound for n1 if n2 = 0
        end
        #####
        for n1=n1bound:n1max # Loop over the index n2
            count += 2 # At this stage (n1,n2) and (-n1,-n2) are two allowed resonances
        end
    end
    #####
    return count # Returning the number of resonance vectors
end


"""
Container of the (n1,n2) resonance vectors to consider

Function that fills in the array of resonance vectors (n1,n2)

ATTENTION, the (n1,n2) are determined for l=lmax
@IMPROVE it would be best to use the same code as in get_nbResVec()
"""
function maketabResVec(lmax::Int64,n1max::Int64)

    # calculate the number
    nbResVec = get_nbResVec(lmax,n1max)
    
    tabResVec = zeros(Int64,2,nbResVec)
    count = 1 # Initialisation of the counter
    #####
    for n2=lmax:-2:0 # Loop over the index n2. ATTENTION, we go by step of 2, and only consider positive values
        #####
        n1bound = -n1max # Value of the bound for n1 if n2 > 0
        #####
        if (n2 == 0)
            n1bound = 1 # Value of the bound for n1 if n2 = 0
        end
        #####
        for n1=n1bound:n1max # Loop over the index n2
            tabResVec[1,count], tabResVec[2,count] = n1, n2 # Adding the resonance (n1,n2)
            #
            count += 1 # Updating the counter
            #####
            tabResVec[1,count], tabResVec[2,count] = -n1, -n2 # Adding the resonance (-n1,-n2)
            #
            count += 1 # Updating the counter
        end
    end
    return tabResVec
end
