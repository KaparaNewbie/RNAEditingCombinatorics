

emptyresults() = DataFrame(
    Fraction=Vector{Float64}(undef, 0),
    FractionRepetition=Vector{Int64}(undef, 0),
    Algorithm=Vector{String}(undef, 0),
    AlgorithmRepetition=Vector{Int64}(undef, 0),
    NumUniqueSamples=Vector{Int64}(undef, 0),
    UniqueSamples=Vector{Vector{String}}(undef, 0),
)



"""
    distinctascending(G)

`G = (V, E)` is a simple graph represented by a dictionary of neighborhood lists.   
(So technically, `V = keys(G)` and `E = {(u, v) | v ∈ G[u]}`.)  
First, sort the vertices in ascending order of degree (within each degree, sort the vertices at random).  
Then, initialize a subset of vertices `V2` and, iteratively, add vertices from `V` to `V2` if they aren't neighbors
of any other vertices already present in `V2` (according to the original edges `E`).  
Return the vertices in `V2` in a vector. Use `sortresults=true` to return them sorted in ascending order of vertices' names.
"""
# @timeit to function distinctascending(G, sortresults::Bool=false)
function distinctascending(G, sortresults::Bool=false)

    # create a vector with tuples of (vertex, degree)
    verticesdegs = [(v, length(G[v])) for v ∈ keys(G)]
    # sort the vector in ascending order of degree; within each degree, sort the vertices at random
    ThreadsX.sort!(
        verticesdegs,
        lt=(x, y) -> (x[2] < y[2] || (x[2] == y[2] && rand([0, 1]) == 1)),
    )

    # iteratively, add vertices from V to V2 - 
    # but only if they aren't neighbors of any other vertices already present in V2 according to E
    V2 = Set{eltype(keys(G))}()
    for (u, deg) ∈ verticesdegs
        if deg == 0 || nofriendships(G, u, V2)
            push!(V2, u)
        end
    end

    # return V2 as a vector (optionally sorted in ascending order of names)
    if sortresults
        return ThreadsX.sort([x for x ∈ V2])
    else
        return collect(V2)
    end

end




"""
    ifriendships(G, u, V2)

Iterate over the possible friendships of `u` and any `v ∈ V2 ⊂ V` according to `G = (V, E)`.  
Technically speaking, `G` is a simple graph represented by a dictionary of neighborhood lists.  
"""
ifriendships(G, u, V2) = imap(v -> u ∈ G[v], V2)


"""
    nofriendships(G, u, V2)

Return `true` if `u` isn't a neighbor of any `v ∈ V2 ⊂ V` according to `G = (V, E)`.  
Technically speaking, `G` is a simple graph represented by a dictionary of neighborhood lists.  
"""
function nofriendships(G, u, V2)
    for friendship ∈ ifriendships(G, u, V2)
        friendship == true && return false
    end
    return true
end




"""
    distinctdescending(G)

`G = (V, E)` is a simple graph represented by a dictionary of neighborhood lists.   
(So technically, `V = keys(G)` and `E = {(u, v) | v ∈ G[u]}`.)  
Iteratively, remove from the graph a random highest-degree vertex and its edges until no edges remain.   
Return the remaining vertices in a vector. Use `sort=true` to return them in ascending order of vertices' names. 

Same as `distinctdescending_old` except it's supposed to be faster due to better management of `maxdeg` & `maxdegvertices`.
"""
# @timeit to function distinctdescending(G, sortresults::Bool=false)
function distinctdescending(G, sortresults::Bool=false)

    G = deepcopy(G) # todo should we remove the deepcopy? it depends on the repetion scheme

    descdegs = sort(unique(length.(values(G))), rev=true) # unique degrees of G, sorted by descending order
    # descdegs = sort(unique(length(vs) for vs ∈ values(G)), rev=true)

    degseltype = typeof(descdegs[1])
    verticeseltype = eltype(keys(G))

    verticesbydeg = Dict{degseltype,Set{verticeseltype}}()

    for (u, vs) ∈ G
        deg = length(vs)
        if !haskey(verticesbydeg, deg)
            verticesbydeg[deg] = Set{verticeseltype}()
        end
        push!(verticesbydeg[deg], u)
    end

    maxdeg = descdegs[1]

    while maxdeg > 0

        u = rand(verticesbydeg[maxdeg])

        # `u` is the single vertex in `verticesbydeg[maxdeg]`
        if length(verticesbydeg[maxdeg]) == 1
            # we can erase the `verticesbydeg[maxdeg]` set altogether
            delete!(verticesbydeg, maxdeg)
            # we also remove the first element from `descdegs`, which is `maxdeg`
            deleteat!(descdegs, 1)
            # there's at least one more vertex other than `u` in `verticesbydeg[maxdeg]`,
        else
            # we only remove `u` from `verticesbydeg[maxdeg]`
            delete!(verticesbydeg[maxdeg], u)
        end

        # update the graph & other variables that store information about degrees

        #= 
        todo 
        when iterating over all `v ∈ G[u]`, instead of updating `verticesbydeg` & `descdegs` separately 
        for each `v` according to its degrees `olddeg` & `newdeg`, we could only track the needed changes, and
        resolve them together after that.
        for example, consider the following case for `v1` and `v2`:

            before:

            olddegv1 = 8
            olddegv2 = 7
            verticesbydeg_A = Dict(8 => {v1}, 7 => {v2}) 

            after:

            newdegv1 = 7
            newdegv2 = 6
            verticesbydeg_C = Dict(7 => {v1}, 6 => {v2})

        if we meet v2 before v1, we first do

            verticesbydeg_B = Dict(8 => {v1}, 6 => {v2}) 

        and only then 

            verticesbydeg_C = Dict(7 => {v1}, 6 => {v2})

        but by knowing all the movements in advance, we could (maybe, for some updates resulting from rand `u` selection):
        1. save some sets creations in `verticesbydeg`, 
        2. and, more importantly, avoid deletions and insertions in `descdegs`
        =#

        for v ∈ G[u]

            olddeg = length(G[v])
            newdeg = olddeg - 1

            # find location of `olddeg` in `descdegs`
            olddegidx = searchsortedlast(descdegs, olddeg, rev=true)

            # remove the edge (u, v) from the graph
            delete!(G[v], u)

            # `v` is the single vertex in `verticesbydeg[olddeg]`
            if length(verticesbydeg[olddeg]) == 1
                # we can erase the `verticesbydeg[olddeg]` set altogether
                delete!(verticesbydeg, olddeg)
                # there's at least one more vertex other than `v` whose degree is `newdeg` (`v`'s updated degree)
                if haskey(verticesbydeg, newdeg)
                    # add `v` to those vertices whose degree is `newdeg`
                    push!(verticesbydeg[newdeg], v)
                    # remove the element at index `olddegidx` from `descdegs`, which is `olddeg`
                    deleteat!(descdegs, olddegidx)
                    # after removing the edge (u, v), `v` is the single vertex whose degree is `newdeg`
                else
                    # create a new set of vertices with degree `newdeg`, and `v` as its single vertex
                    verticesbydeg[newdeg] = Set{verticeseltype}()
                    push!(verticesbydeg[newdeg], v)
                    # as `v` was the only vertex with `olddeg` before removing (u, v),
                    # and is the only vertex with degree `newdeg` after that,
                    # we simply replace `olddeg` with `newdeg` in `descdegs`
                    descdegs[olddegidx] = newdeg
                end
                # there remain other vertices whose degree is `olddeg`
            else
                # remove `v` from `verticesbydeg[olddeg]`
                delete!(verticesbydeg[olddeg], v)
                # there's at least one more vertex other than `v` whose degree is `newdeg` (`v`'s updated degree)
                if haskey(verticesbydeg, newdeg)
                    # add `v` to those vertices whose degree is `newdeg`
                    push!(verticesbydeg[newdeg], v)
                    # after removing the edge (u, v), `v` is the single vertex whose degree is `newdeg`
                else
                    # create a new set of vertices with degree `newdeg`, and `v` as its single vertex
                    verticesbydeg[newdeg] = Set{verticeseltype}()
                    push!(verticesbydeg[newdeg], v)
                    # as there remain other vertices whose degree is `olddeg` after removing (u, v),
                    # but `v` is the only vertex with degree `newdeg` after that,
                    # we add `newdeg` to `descdegs` at index `olddegidx + 1`
                    insert!(descdegs, olddegidx + 1, newdeg)
                end
            end

        end

        # remove `u` from the graph
        delete!(G, u)

        # update the maximum degree in the updated graph
        maxdeg = descdegs[1]

    end

    # return the remaining vertices as a vector, optionally sorted in ascending order of names
    if sortresults
        return ThreadsX.sort([x for x ∈ keys(G)])
    else
        return [x for x ∈ keys(G)]
    end

end


function distinctunordered(G, sortresults::Bool=false)
    # iterate over a random permutation of `V` and, 
    # for each `u` in that permutation, add it to `V2` -
    # but only if `u` isn't a neighbor of any of the other vertices already present in V2 
    # (according to E)
    V = shuffle(collect(keys(G)))
    V2 = Set{eltype(keys(G))}()
    for u ∈ V
        if nofriendships(G, u, V2)
            push!(V2, u)
        end
    end

    # return V2 as a vector (optionally sorted in ascending order of names)
    if sortresults
        return ThreadsX.sort([x for x ∈ V2])
    else
        return collect(V2)
    end

end


const algfuncs = Dict(
    "Ascending" => distinctascending,
    "Descending" => distinctdescending,
    "Unordered" => distinctunordered
)


"""
Find subsets of distinct unique samples (unique reads/proteins).
"""
function solve(
    sampleG::Dict,
    fraction::Float64,
    fracrepetition::Int64,
    algrepetitions::Int64,
    threaded::Bool=false,
    sortresults::Bool=false,
    algs::Vector{String}=["Ascending", "Descending", "Unordered"],
)
    @info "$(loggingtime())\tsolve" fraction fracrepetition algrepetitions threaded algs algfuncs sortresults myid()
    # run the different algorithms
    mapinputs = [
        (alg, algrepetition)
        for alg ∈ algs
        for algrepetition ∈ 1:algrepetitions
    ]
    mapfunc = threaded ? ThreadsX.map : map
    # @timeit to "algfunc (obtaining distinct isoforoms using ascending / aescending function)" begin
    #     distinctsamples_alg_algrepetition_sets = mapfunc(mapinputs) do (alg, algrepetition)
    #         algfunc = algfuncs[alg] # the function implementing `alg`
    #         distinctsamples = algfunc(sampleG, sortresults)  # distinct unique samples obtained by `algfunc`
    #         (distinctsamples, alg, algrepetition)
    #     end    
    # end
    distinctsamples_alg_algrepetition_sets = mapfunc(mapinputs) do (alg, algrepetition)
        algfunc = algfuncs[alg] # the function implementing `alg`
        distinctsamples = algfunc(sampleG, sortresults)  # distinct unique samples obtained by `algfunc`
        (distinctsamples, alg, algrepetition)
    end


    # create an empty results df 
    results = emptyresults()
    # add the distinct unique samples (reads/proteins) to the results
    for (distinctsamples, alg, algrepetition) ∈ distinctsamples_alg_algrepetition_sets
        push!(
            results,
            [fraction, fracrepetition, alg, algrepetition, length(distinctsamples), distinctsamples]
        )
    end

    # concatanate the vector of distinct unique samples in each result into a semicolon-separated string
    transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)
    # return results
    return results
end

