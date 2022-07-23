

"""
# We define a list of the (np,nq) for which the response matrix must be computed
# It is useful to :
# + Use the symmetry (np,nq)
# + Efficiently parallelise on large core numbers
# Function that initialises the table of (np,nq)
# for which the response matrix must be computed
"""
function makeTabnpnq(nradial::Int64)
    nb_npnq = convert(Int64,nradial*(nradial+1)/2)
    tab_npnq = zeros(Int64,2,nb_npnq)
    count = 1 # Counter of the position in the array
    #####
    for np=1:nradial # Loop over the index np
        for nq=np:nradial # Loop over the index nq. ATTENTION, we use the symmetry (np,nq) <-> (nq,np) to reduce the number of elements to consider
            tab_npnq[1,count], tab_npnq[2,count] = np,nq # Filling in the array
            count += 1 # Updating the counter
        end
    end
    return tab_npnq
end
