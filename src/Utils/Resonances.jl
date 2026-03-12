"""
helper functions for making resonance containers

example:
lharmonic=2
n1max=2
nbResVec = get_nbResVec(lharmonic,n1max) # Number of resonance vectors. ATTENTION, it is for the harmonics lharmonic
tabResVec = maketabResVec(lharmonic,n1max) # Filling in the array of resonance vectors (n1,n2)
println(tabResVec)

"""


function GetNbResVec(lharmonic::Int64,n1max::Int64,dimension::Int64=3,isVanishing::Bool=true)
    if dimension == 2
        return GetNbResVec2d(lharmonic,n1max,isVanishing)
    elseif dimension == 3
        return GetNbResVec3d(lharmonic,n1max)
    else
        error("Unknow dimension in GetNbResVec")
    end
end

function MakeTabResVec(lharmonic::Int64,n1max::Int64,dimension::Int64=3,isVanishing::Bool=true)
    if dimension == 2
        return MakeTabResVec2d(lharmonic,n1max,isVanishing)
    elseif dimension == 3
        return MakeTabResVec3d(lharmonic,n1max)
    else
        error("Unknow dimension in MakeTabResVec")
    end
end


"""
    GetNbResVec3d(lharmonic::Int64, n1max::Int64)

Calculate the total number of resonances (bn = (n1, n2)) to consider for the harmonics `lharmonic`.

This function calculates the number of resonance vectors (bn) that satisfy the following constraints:
- |n2| <= lharmonic
- (lharmonic - n2) must be even
- |n1| <= n1max
- (n1, n2) = (0, 0) does not contribute

# Arguments:
- `lharmonic::Int64`: Harmonics for which the resonances are calculated.
- `n1max::Int64`: Maximum value for the component n1.

# Returns:
- `Int64`: Total number of allowed resonance vectors.

# Example:
```julia
lharmonic = 3
n1max = 2
resonance_count = GetNbResVec3d(lharmonic, n1max)  # Returns the number of resonances for lharmonic=3 and n1max=2
```
"""
function GetNbResVec3d(lharmonic::Int64,n1max::Int64)
    count = 0 # Initialisation of the counter
    #####
    for n2=lharmonic:-2:0 # Loop over the index n2. ATTENTION, we go by step of 2, and only consider positive values
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

ATTENTION, the (n1,n2) are determined for l=lharmonic
@IMPROVE it would be best to use the same code as in get_nbResVec()
"""
function MakeTabResVec3d(lharmonic::Int64,n1max::Int64)

    # calculate the number
    nbResVec = GetNbResVec(lharmonic,n1max,3)

    tabResVec = zeros(Int64,2,nbResVec)
    count = 1 # Initialisation of the counter
    #####
    for n2=lharmonic:-2:0 # Loop over the index n2. ATTENTION, we go by step of 2, and only consider positive values
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
    return nbResVec, tabResVec
end


"""
# Function that returns the total number
# of resonances bn=(n1,n2) to consider
# for the harmonics lharmonic for discs, 
# assuming a L-integration over [0, inf).
# There a few constraints to satisfy:
# + n2 = ±lharmonic for each value of n1 if isVanishing=false
# + n2 = +lharmonic for each value of n1 if isVanishing=true
# + |n1| <= n1max
# + (n1,n2) = (0,0) does not contribute
# ATTENTION, the (n1,n2) are determined for l=lharmonic
"""
function GetNbResVec2d(lharmonic::Int64,n1max::Int64,isVanishing::Bool=true)
    count = 0 # Initialisation of the counter
    #####
    # In the the 2D, this package also integrates L over [0, inf).
    # However, the response matrix should integrate L over (-inf, inf), since the L variable in the 2D case stands for Lz=R*vphi

    # If the DF is even or odd in L, then :
    # We can keep the integration of L over [0, inf), but at the cost of adding the contribution of n2=-lharmonic to the resonance summation, for each values of n1

    # If the DF vanishes smoothly at L=0, then:
    # We can keep the integration of L over [0, inf), with only a contribbution from n2=+lharmonic to the resonance summation, for each values of n1

    #####
    if lharmonic != 0
        for n1=-n1max:n1max # Loop over the index n1
            if (isVanishing)
                count += 1 # Updating the counter : For each n1, we have the contributions from n2=lharmonic
            else
                count += 2 # Updating the counter : For each n1, we have the contributions from n2=-lharmonic and n2=lharmonic
            end
        end
    else
        for n1=1:n1max # Loop over the index n1
            if (isVanishing)
                count += 2 # Updating the counter (n1,lharmonic), and (-n1,lharmonic)
            else
                count += 4 # Updating the counter (n1,lharmonic), (n1,-lharmonic), (-n1,lharmonic) and (-n1,-lharmonic)   
            end
        end
    end
    #####
    return count # Returning the number of resonance vectors
end


"""
Container of the (n1,n2) resonance vectors to consider

Function that fills in the array of resonance vectors (n1,n2)

ATTENTION, the (n1,n2) are determined for l=lharmonic
@IMPROVE it would be best to use the same code as in get_nbResVec()
"""
function MakeTabResVec2d(lharmonic::Int64,n1max::Int64,isVanishing::Bool=true)
    # calculate the number
    nbResVec = GetNbResVec(lharmonic,n1max,2,isVanishing)

    tabResVec = zeros(Int64,2,nbResVec)
    count = 1 # Initialisation of the counter
    #####
    # In the the 2D, this package also integrates L over [0, inf).
    # However, the response matrix should integrate L over (-inf, inf), since the L variable in the 2D case stands for Lz=R*vphi

    # If the DF is even or odd in L, then :
    # We can keep the integration of L over [0, inf), but at the cost of adding the contribution of n2=-lharmonic to the resonance summation, for each values of n1

    # If the DF vanishes smoothly at L=0, then:
    # We can keep the integration of L over [0, inf), with only a contribbution from n2=+lharmonic to the resonance summation, for each values of n1

    #####
    if lharmonic != 0
        for n1=-n1max:n1max # Loop over the index n1
            tabResVec[1,count], tabResVec[2,count] = n1, lharmonic # Adding the resonance (n1,lharmonic)
            #
            count += 1 # Updating the counter
            #
            if (!isVanishing)
                tabResVec[1,count], tabResVec[2,count] = n1, -lharmonic # Adding the resonance (n1,-lharmonic)
                #
                count += 1 # Updating the counter
            end
        end
    else
        for n1=1:n1max # Loop over the index n1
            tabResVec[1,count], tabResVec[2,count] = n1, lharmonic # Adding the resonance (n1,lharmonic)
            #
            count += 1 # Updating the counter
            #
            tabResVec[1,count], tabResVec[2,count] = -n1, lharmonic # Adding the resonance (-n1,lharmonic)
            #
            count += 1 # Updating the counter
            #
            if (!isVanishing)
                tabResVec[1,count], tabResVec[2,count] = n1, -lharmonic # Adding the resonance (n1,-lharmonic)
                #
                count += 1 # Updating the counter
                #
                tabResVec[1,count], tabResVec[2,count] = -n1, -lharmonic # Adding the resonance (-n1,-lharmonic)
                #
                count += 1 # Updating the counter
            end
        end
    end

    # For the even DF case for lmax=0, the resonance number appear twice, which is accurate
    # Indeed, lmax=0 implied n2=0, and in that case the integration over (-inf, inf) is simply twice that over (0, inf)
    # We might be able to optimize this, but the current form allows us not to change anything in the package 

    return nbResVec, tabResVec
end
