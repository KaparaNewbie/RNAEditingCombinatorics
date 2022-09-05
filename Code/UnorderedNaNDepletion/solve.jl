



# function solve(
#     sampleG::Dict,
#     fraction::Float64,
#     fracrepetition::Int64,
#     algrepetitions::Int64
# )
#     @info "$(loggingtime())\tsolve" fraction fracrepetition algrepetitions myid()

#     # find subsets of distinct unique samples (unique reads/proteins),using ascending & descending algorithms 
#     # in $repetions in parallel

#     # todo define a parameter that determinnes if both algs are spawned and then fetched, 
#     # or if the second is spawned only after the first is finished (the way it's implemented now)

#     distinctasc = [@spawn distinctascending(sampleG) for _ ∈ 1:algrepetitions]
#     distinctasc = [fetch(t) for t ∈ distinctasc] # unique samples obtained by the ascending algorithm
#     distinctdesc = [@spawn distinctdescending(sampleG) for _ ∈ 1:algrepetitions]
#     distinctdesc = [fetch(t) for t ∈ distinctdesc] # unique samples obtained by the descending algorithm

#     # add the compatible unique samples to the results
#     results = emptyresults()
#     for (algrepetition, distinctsamples) ∈ enumerate(distinctasc)
#         push!(results, [fraction, fracrepetition, "Ascending", algrepetition, length(distinctsamples), distinctsamples])
#     end
#     for (algrepetition, distinctsamples) ∈ enumerate(distinctdesc)
#         push!(results, [fraction, fracrepetition, "Descending", algrepetition, length(distinctsamples), distinctsamples])
#     end
#     # concatanate the vector of distinct unique samples in each result into a semicolon-separated string
#     transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)
#     # return results
#     return results
# end



# "Sequential version of solve. Each alg is run only once."
# function solve(
#     sampleG::Dict,
#     fraction::Float64,
#     fracrepetition::Int64
# )
#     @info "$(loggingtime())\tsolve" fraction fracrepetition myid()
#     # find subsets of distinct unique samples (unique reads/proteins) using ascending & descending algorithms 
#     distinctasc = distinctascending(sampleG) # distinct unique samples obtained by the ascending algorithm
#     distinctdesc = distinctdescending(sampleG) # distinct unique samples obtained by the descending algorithm
#     # create an empty results df 
#     results = emptyresults()
#     # add the distinct unique samples to the results
#     algrepetition = 1
#     push!(results, [fraction, fracrepetition, "Ascending", algrepetition, length(distinctasc), distinctasc])
#     push!(results, [fraction, fracrepetition, "Descending", algrepetition, length(distinctdesc), distinctdesc])
#     # concatanate the vector of distinct unique samples in each result into a semicolon-separated string
#     transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)
#     # return results
#     return results
# end



# struct SolveThreaded{x} end
# SolveThreaded(x) = SolveThreaded{x}()



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
    # algfuncs::Dict{String,Function}=Dict(
    #     "Ascending" => distinctascending,
    #     "Descending" => distinctdescending,
    #     "Unordered" => distinctunordered
    # )
)
    @info "$(loggingtime())\tsolve" fraction fracrepetition algrepetitions threaded algs algfuncs sortresults myid()
    # run the different algorithms
    mapinputs = [
        (alg, algrepetition)
        for alg ∈ algs
        for algrepetition ∈ 1:algrepetitions
    ]
    mapfunc = threaded ? ThreadsX.map : map
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










"""
    nofriendships_2(G, u, V2)

Return `true` if `u` isn't a neighbor of any `v ∈ V2 ⊂ V` according to `G = (V, E)`.  
Technically speaking, `G` is a simple graph represented by a dictionary of neighborhood lists.  
This is a threaded version of `nofriendships`.
"""
# nofriendships_2(G, u, V2) = ThreadsX.findfirst([u ∈ G[v] for v ∈ V2]) === nothing
nofriendships_2(G, u, V2) = !ThreadsX.any([u ∈ G[v] for v ∈ V2])



function distinctunordered_2(G, sortresults::Bool=false)
    # iterate over a random permutation of `V` and, 
    # for each `u` in that permutation, add it to `V2` -
    # but only if `u` isn't a neighbor of any of the other vertices already present in V2 
    # (according to E)
    V = shuffle(collect(keys(G)))
    degrees = [length(G[u]) for u ∈ V]
    V2 = Set{eltype(V)}()
    for (u, deg) ∈ zip(V, degrees)
        if deg == 0 || nofriendships_2(G, u, V2)
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



function distinctascending_2(G, sortresults::Bool=false)

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
        if deg == 0 || nofriendships_2(G, u, V2)
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
    distinctdescending_old(G)

`G = (V, E)` is a simple graph represented by a dictionary of neighborhood lists.   
(So technically, `V = keys(G)` and `E = {(u, v) | v ∈ G[u]}`.)  
Iteratively, remove from the graph a random highest-degree vertex and its edges until no edges remain.   
Return the remaining vertices in a vector. Use `sortresults=true` to return them in ascending order of vertices' names. 
"""
function distinctdescending_old(G, sortresults::Bool=false)

    G = deepcopy(G) # todo should we remove the deepcopy? it depends on the repetion scheme

    degrees = Dict(v => length(G[v]) for v ∈ keys(G))
    maxdeg = ThreadsX.maximum(values(degrees)) # the highest degree in the original graph
    maxdegvertices = [v for (v, deg) ∈ degrees if deg == maxdeg]

    # remove a random highest-degree vertex and its edges until no edges remain
    # maxdeg = 1 # artificial value in order to run the loop at least once
    while maxdeg > 0
        # find a random highest-degree vertex u
        # maxdeg = ThreadsX.maximum(values(degrees)) # the highest degree in the original graph
        # maxdegvertices = [v for (v, deg) ∈ degrees if deg == maxdeg]
        u = rand(maxdegvertices)
        # update the graph by removing (u, v) edges and the random highest-degree vertex u itself
        for v ∈ G[u]  # todo parallelize this loop? `degrees` will have to be a ThreadSafeDict https://github.com/wherrera10/ThreadSafeDicts.jl
            delete!(G[v], u)
            degrees[v] -= 1
        end
        delete!(G, u)
        delete!(degrees, u)
        # update the maximum degree in the updated graph
        maxdeg = ThreadsX.maximum(values(degrees))
        maxdegvertices = [v for (v, deg) ∈ degrees if deg == maxdeg]
    end

    # return the remaining vertices as a vector, optionally sorted in ascending order of names
    if sortresults
        return ThreadsX.sort([x for x ∈ keys(G)])
    else
        return collect(keys(G))
    end
end




# infile = "D.pealeii/MpileupAndTranscripts/IlluminaExample2/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv"
# samplename = "comp141044_c0_seq2"
# idcol = "Protein"
# outdir = "D.pealeii/MpileupAndTranscripts/IlluminaExample2"
# fracstep = 1.0
# maxfrac = 1.0
# fracrepetitions = 1
# algrepetitions = 1
# delim = "\t"
# datatype = "Proteins"
# fraction = 1.0
# fracrepetition = 1
# algrepetitions = 10
# firstcolpos = 15
# testfraction = 0.06

# df, firstcolpos = preparedf!(infile, delim, datatype, idcol, firstcolpos, testfraction)

# G = indistinguishable_rows(df, idcol; firstcolpos)

# # @benchmark distinctunordered($G)
# # @benchmark distinctascending($G)
# # @benchmark distinctdescending($G)

# # @benchmark distinctunordered_2($G) 
# # @benchmark distinctascending_2($G) 


# fraction = 1.0
# nrows = size(df, 1)
# nsamplerows = convert(Int, round(fraction * nrows)) 

# # sampleG = get_graph_sample(G, nsamplerows, df, idcol, 1892)

# # fracrepetition = 1 
# # algrepetitions = 5
# # threaded = true
# # results = solve(
# #     sampleG, fraction, fracrepetition, algrepetitions, threaded
# # )

# allresults = map(1:5) do i
#     sampleG = get_graph_sample(G, nsamplerows, df, idcol, i)
#     fracrepetition = i
#     algrepetitions = 10
#     threaded = true
#     solve(sampleG, fraction, fracrepetition, algrepetitions, threaded)
# end
# results = vcat(allresults...)
# sort!(results, "Algorithm")


# jaccardindex(s1, s2) = length(s1 ∩ s2) / length(s1 ∪ s2)

# function jaccrdmatrix(sets)
#     jm = Matrix{Float64}(undef, length(sets), length(sets))
#     for x ∈ eachindex(sets)
#         s1 = sets[x]
#         for y ∈ eachindex(sets)
#             s2 = sets[y]
#             jm[x, y] = jaccardindex(s1, s2)
#         end
#     end
#     jm
# end


# using CairoMakie


# sets = Set.(split.(results[:, "UniqueSamples"], ","))
# jm = jaccrdmatrix(sets)

# # ticks = select(results, "Algorithm" => ByRow(x -> x[1:1]))[:, 1]
# ticks = results[:, "Algorithm"]
# l = length(ticks) 
# mids = [25, 75, 125]
# ticks = map(eachindex(ticks)) do i
#     i ∈ mids ? ticks[i] : ""
# end
# xticks = ticks;
# yticks = ticks;

# fig = Figure(resolution = (1200, 1200), fontsize = 20)
# ax = Axis(fig[1, 1], xticks = (1:l, xticks), yticks = (1:l, yticks))
# hmap = heatmap!(ax, jm, colormap = :plasma)
# Colorbar(fig[1, 2], hmap; label = "Jaccard index", width = 25, ticksize = 15, labelsize=20)
# # ax.xticklabelrotation = π / 3
# # ax.xticklabelalign = (:right, :center)
# display(fig)



# using MultivariateStats, RDatasets, Plots

# # load iris dataset
# iris = dataset("datasets", "iris")

# # Performing PCA on Iris data set:

# Xtr = jm[1:2:end,:]'
# Xtr_labels = results[1:2:end, "Algorithm"]
# # split other half to testing set
# Xte = jm[2:2:end,:]'
# Xte_labels = results[2:2:end, "Algorithm"]

# # Suppose Xtr and Xte are training and testing data matrix, with each observation in a column. 
# # We train a PCA model, allowing up to 3 dimensions:

# M = fit(PCA, Xtr; maxoutdim=3)

# # Then, apply PCA model to the testing set

# Yte = predict(M, Xte)

# # And, reconstruct testing observations (approximately) to the original space

# Xr = reconstruct(M, Yte)  # does this line affect the plot?

# # Now, we group results by testing set labels for color coding and visualize first 3 principal components in 3D plot

# asc = Yte[:,Xte_labels.=="Ascending"]
# desc = Yte[:,Xte_labels.=="Descending"]
# unord = Yte[:,Xte_labels.=="Unordered"]

# p = Plots.scatter(asc[1,:],asc[2,:],asc[3,:],marker=:circle,linewidth=0,label="Ascending")
# Plots.scatter!(desc[1,:],desc[2,:],desc[3,:],marker=:circle,linewidth=0,label="Descending")
# Plots.scatter!(unord[1,:],unord[2,:],unord[3,:],marker=:circle,linewidth=0,label="Unordered")
# Plots.plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")

