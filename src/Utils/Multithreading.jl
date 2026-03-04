# Suppose we contruct a vector of local quantities for multi-threading
# e.g. list_tab = [tab for i=1:Threads.nthreads()]
# get_pos_threadid() returns the dictionary : threadid => index
# Such that the following code works
# pos_threadid = get_pos_threadid()
# Threads.@threads :static for i=1:n
#     tid = Threads.threadid()
#     id = pos_threadid(tid)
#     foo!(list_tab[id])
# end
function get_pos_threadid()

    list_threadid = zeros(Int64, Threads.nthreads())

    Threads.@threads :static for i=1:Threads.nthreads()
        list_threadid[i] = Threads.threadid()
    end

    pos_threadid = Dict(tid => i for (i, tid) in enumerate(list_threadid))

    return pos_threadid

end