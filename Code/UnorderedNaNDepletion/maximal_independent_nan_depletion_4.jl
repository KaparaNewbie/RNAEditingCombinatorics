# this is an old version of maximal_independet_set_5.jl  
# algorithmically speaking, everything in here exists in the new file - 
# this version is only kept because some runs were made with it


# https://discourse.julialang.org/t/julia-python-equivalent-of-main/35433
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
else
    using BenchmarkTools
end

using ArgParse
using CSV
using DelimitedFiles
using DataFrames
using LinearAlgebra
using StatsBase  # for StatsBase.sample
using IterTools  # for IterTools.imap
import Base.Threads.@spawn
using Random
using BioSequences
using Transducers
# using Distributed
# using Dagger


const ∅ = Set()


"""
    indistinguishable_rows(M, ids)

Construct a simple graph `G` represented by a dictionary of neighborhood lists.  
Each row of `M` (a unique read) is a vertex in the graph.  
Two unique reads `i` and `j`, `i ≠ j`, are connected (`j ∈ G[i]`) if they cannot be distinguished by any editing position 
(the first of which is the column whose index is `firstcolpos`). 
That is, `i` and `j` are uncompatible and thus cannot be considered both as unique reads in final analysis. 
Else, `j ∉ G[i]`.
"""
function indistinguishable_rows(M::Matrix{Set{AminoAcid}}, ids)
    # informativecols = tcollect(i for (i, col) in enumerate(eachcol(M)) if length(unique(col)) > 1)
    # M = M[:, informativecols]
    nrows = size(M, 1)
    length(ids) == length(unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )
    nodeseltype = eltype(ids)
    AAsets = Set([x for row ∈ eachrow(M) for x ∈ row])
    emptyAAintersections = Set(
        [(x, y)
         for x ∈ AAsets for y ∈ AAsets
         if x ∩ y == ∅]
    )
    Gs = tcollect(indistinguishable_vs_for_u(M, emptyAAintersections, ids, nodeseltype, i) for i ∈ 1:nrows)
    # G = merge(Gs...)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # each g has a single (k, v) pair
    return G
end


function indistinguishable_vs_for_u(
    M::Matrix{Set{AminoAcid}}, emptyAAintersections, ids, nodeseltype, i
)
    u = ids[i]
    @views rowᵢ = M[i, :]
    G = Dict(u => Set{nodeseltype}())
    nrows = size(M, 1)
    for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
        i == j && continue
        v = ids[j]
        @views rowⱼ = M[j, :]
        distinct = false
        for (x, y) ∈ zip(rowᵢ, rowⱼ)
            # if we can find at least one non-missing col for which both rows have different values,
            # then we can consider these two rows as different (distinct)
            if (x, y) ∈ emptyAAintersections
                distinct = true
                break
            end
        end
        # if we can't differ the rows, we consider them as indistinguishable one from each other, 
        # and thus add them to the neighborhood graph
        !distinct && push!(G[u], v)
    end
    return G
end


function indistinguishable_rows(M::Union{Matrix{Float64},Matrix{Int}}, ids)
    # informativecols = tcollect(i for (i, col) in enumerate(eachcol(M)) if length(unique(col)) > 1)
    # M = M[:, informativecols]
    nrows = size(M, 1)
    length(ids) == length(unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )
    nodeseltype = eltype(ids)
    @views M = replace(M, -1 => missing)
    Gs = tcollect(indistinguishable_vs_for_u(M, ids, nodeseltype, i) for i ∈ 1:nrows)
    # G = merge(Gs...)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g)
    return G
end


function indistinguishable_vs_for_u(
    M::Union{Matrix{Union{Missing,Int64}},Matrix{Union{Missing,Float64}}}, ids, nodeseltype, i
)
    u = ids[i]
    @views rowᵢ = M[i, :]
    G = Dict(u => Set{nodeseltype}())
    nrows = size(M, 1)
    for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
        i == j && continue
        v = ids[j]
        @views rowⱼ = M[j, :]
        distinct = false
        for (x, y) ∈ zip(rowᵢ, rowⱼ)
            # if we can find at least one non-missing col for which both rows have different values,
            # then we can consider these two rows as different (distinct)
            if x !== missing && y !== missing && x != y
                distinct = true
                break
            end
        end
        # if we can't differ the rows, we consider them as indistinguishable one from each other, 
        # and thus add them to the neighborhood graph
        !distinct && push!(G[u], v)
    end
    return G
end



function indistinguishable_rows(df::DataFrame, idcol; firstcolpos::Int=2, areuniuqe::Bool=false)
    if !areuniuqe
        df = unique(df, idcol)
        areuniuqe = true
    end
    M = Matrix(df[:, firstcolpos:end])
    ids = df[:, idcol]
    return indistinguishable_rows(M, ids)
end


"""
    uncompatiblerows(M, firstcolpos)

Construct a simple graph `G` represented by a dictionary of neighborhood lists.  
Each row of `M` (a unique read) is a vertex in the graph.  
Two unique reads `i` and `j`, `i ≠ j`, are connected (`j ∈ G[i]`) if they cannot be distinguished by any editing position 
(the first of which is the column whose index is `firstcolpos`). 
That is, `i` and `j` are uncompatible and thus cannot be considered both as unique reads in final analysis. 
Else, `j ∉ G[i]`.
"""
# function uncompatiblerows(M::Union{Matrix{Float64},Matrix{Int}}, ids; areuniuqe::Bool=false)
#     if !areuniuqe
#         @views M = unique(M, dims=1)
#         ids = unique(ids)
#     end
#     nrows = size(M, 1)
#     @assert length(ids) == nrows
#     idseltype = eltype(ids)
#     G = Dict(u => Set{idseltype}() for u ∈ ids)
#     @views M = replace(M, -1 => missing)
#     Threads.@threads for i ∈ 1:nrows
#         @views rowᵢ = M[i, :]
#         u = ids[i]
#         Threads.@threads for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
#             i == j && continuek
#             @views rowⱼ = M[j, :]
#             v = ids[j]
#             compatible = false
#             for (x, y) ∈ zip(rowᵢ, rowⱼ)
#                 # if we can find at least one non-missing col for which both rows have different values,
#                 # then we can consider these two rows as different (compatible)
#                 if x !== missing && y !== missing && x != y
#                     compatible = true
#                     break
#                 end
#             end
#             # if we can't differ the rows, we consider them as uncompatible one to each other, 
#             # and thus add them to the neighborhood graph
#             compatible == false && push!(G[u], v)
#         end
#     end
#     return G
# end


# function uncompatiblerows(M::Matrix{Set{AminoAcid}}, ids; areuniuqe::Bool=false)
#     if !areuniuqe
#         @views M = unique(M, dims=1)
#         ids = unique(ids)
#     end
#     nrows = size(M, 1)
#     @assert length(ids) == nrows
#     idseltype = eltype(ids)
#     G = Dict(u => Set{idseltype}() for u ∈ ids)
#     AAsets = Set([x for row ∈ eachrow(M) for x ∈ row])
#     emptyAAintersections = Set(
#         [(x, y)
#          for x ∈ AAsets for y ∈ AAsets
#          if x ∩ y == ∅]
#     )
#     Threads.@threads for i ∈ 1:nrows
#         @views rowᵢ = M[i, :]
#         u = ids[i]
#         Threads.@threads for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
#             i == j && continue
#             @views rowⱼ = M[j, :]
#             v = ids[j]
#             compatible = false
#             for (x, y) ∈ zip(rowᵢ, rowⱼ)
#                 # if we can find at least one col for which the rows' intersection is empty,
#                 # that is, if the rows don't share any possible amino acid for that position,
#                 # then we can consider these two rows as different (compatible)
#                 if (x, y) ∈ emptyAAintersections
#                     compatible = true
#                     break
#                 end
#             end
#             # if we can't differ the rows, we consider them as uncompatible one to each other, 
#             # and thus add them to the neighborhood graph
#             compatible == false && push!(G[u], v)
#         end
#     end
#     return G
# end


# function uncompatiblerows(df::DataFrame, idcol; firstcolpos::Int=2, areuniuqe::Bool=false)
#     if !areuniuqe
#         df = unique(df, idcol)
#         areuniuqe = true
#     end
#     return uncompatiblerows(Matrix(df[:, firstcolpos:end]), df[:, idcol]; areuniuqe)
# end



"""
    subgraph(G, V2)

Return a subgraph `G2 ⊆ G` induced by the vertices `V2 ⊆ V`, where `V` are the original vertices of `G`.
"""
function subgraph(G, V2)
    V2eltype = eltype(V2)
    G2 = Dict(u => Set{V2eltype}() for u ∈ V2)
    length(G2) == length(V2) || error("The vertices in V2 should be unique.")
    Threads.@threads for u ∈ V2
        Threads.@threads for v ∈ V2 # compare vertices u & v if u != v
            u == v && continue
            if v ∈ G[u]
                push!(G2[u], v)
            end
        end
    end
    return G2
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
    distinctascending(G)

`G = (V, E)` is a simple graph represented by a dictionary of neighborhood lists.   
(So technically, `V = keys(G)` and `E = {(u, v) | v ∈ G[u]}`.)  
First, sort the vertices in ascending order of degree (within each degree, sort the vertices at random).  
Then, initialize a subset of vertices `V2` and, iteratively, add vertices from `V` to `V2` if they aren't neighbors
of any other vertices already present in `V2` (according to the original edges `E`).  
Return the vertices in `V2` in a vector sorted in ascending order of vertices' names.
"""
function distinctascending(G)

    G = deepcopy(G)

    # create a vector with tuples of (vertex, degree)
    verticesdegs = [(v, length(G[v])) for v ∈ keys(G)]
    # sort the vector in ascending order of degree; within each degree, sort the vertices at random
    verticesdegs = sort(
        verticesdegs,
        lt=(x, y) -> (x[2] > y[2] || (x[2] == y[2] && rand([0, 1]) == 1)),
        rev=true
    )

    # iteratively, add vertices from V to V2 - 
    # but only if they aren't neighbors of any other vertices already present in V2 according to E
    # V2 = Set{Integer}()
    V2 = Set{eltype(collect(keys(G)))}()
    for (u, deg) ∈ verticesdegs
        if deg == 0 || nofriendships(G, u, V2)
            push!(V2, u)
        end
    end

    # return V2 as a vector sorted in ascending order of names
    return sort([x for x ∈ V2])

end

"""
    distinctdescending(G)

`G = (V, E)` is a simple graph represented by a dictionary of neighborhood lists.   
(So technically, `V = keys(G)` and `E = {(u, v) | v ∈ G[u]}`.)  
Iteratively, remove from the graph a random highest-degree vertex and its edges until no edges remain.   
Return the remaining vertices in a vector sorted in ascending order of vertices' names.
"""
function distinctdescending(G)

    G = deepcopy(G)

    degrees = Dict(v => length(G[v]) for v ∈ keys(G))
    maxdeg = maximum(values(degrees)) # the highest degree in the original graph

    # remove a random highest-degree vertex and its edges until no edges remain
    while maxdeg > 0
        # find a random highest-degree vertex u
        maxdegvertices = [v for (v, deg) ∈ degrees if deg == maxdeg]
        u = rand(maxdegvertices)
        # update the graph by removing (u, v) edges and the random highest-degree vertex u itself
        for v ∈ G[u]
            delete!(G[v], u)
            degrees[v] -= 1
        end
        delete!(G, u)
        delete!(degrees, u)
        # update the maximum degree in the updated graph
        maxdeg = maximum(values(degrees))
    end

    # return the remaining vertices as a vector sorted in ascending order of names
    return sort([x for x ∈ keys(G)])
end



function parsecmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--infiles"
        help = "One or more csv files representing unique reads."
        nargs = '+'
        action = :store_arg
        required = true
        "--delim"
        help = "Delimiter for input/output csv files."
        arg_type = String
        default = "\t"
        # "--samplesnames"
        # help = "Samples' names of corresponding `infiles`."
        # nargs = '+'
        # action = :store_arg
        # required = true
        "--postfix"
        help = "Postfix of input files. Should be the *full* postfix, e.g. `.unique_reads.csv`."
        required = true
        "--idcol"
        help = "Label of unique samples columns. Typically `Transcripts` for `Reads` and `Protein` for `Proteins`."
        required = true
        "--editingcol"
        help = "Label of column indicating editing level. Any value > 0 indicates the row is edited."
        required = true
        "--firstcolpos"
        help = "Int location of the first editing position column of each file in `infiles`. As of now, should be 9 for `Reads` and 15 for `Proteins`."
        arg_type = Int
        required = true
        "--removeuneditedrows"
        help = "Remove rows that are not edited."
        action = :store_true
        "--allrowsedited"
        help = "Indicate that all rows of each file in `infiles` are edited at least once."
        action = :store_true
        "--datatype"
        help = "Data type of the input files. Either `Reads` or `Proteins`."
        required = true
        "--outdir"
        help = "Write output files to this directory."
        required = true
        "--fracstep"
        help = "Step size for the fraction of the dataset to be used for each iteration, s.t. `fractions = collect(fracstep:fracstep:1.0) ∪ 1.0`."
        arg_type = Float64
        default = 0.1
        "--maxfrac"
        help = "Maximum fraction of the dataset to be sampled."
        arg_type = Float64
        default = 1.0
        "--fracrepetitions"
        help = "Number of repetitions of the sampling procedure for each fraction of the data."
        arg_type = Int
        default = 10
        "--algrepetitions"
        help = "Number of repetitions for each algorithm is run on each repetition of each fraction of the data."
        arg_type = Int
        default = 5
        "--flatparallelism"
        help = "If `true`, don't use nested parallism. Maybe it's actually better?"
        action = :store_true
        "--testfraction"
        help = "Fraction of the dataset to be used for testing."
        arg_type = Union{Float64,Nothing}
        default = nothing
        "--randseed"
        help = "Random seed for sampling test data."
        arg_type = Int
        default = 1892

    end
    return parse_args(s)
end


toAAset(x) = Set(map(aa -> convert(AminoAcid, only(aa)), split(x, ",")))


function preparedf(
    infile, delim,
    removeuneditedrows, allrowsedited, editingcol,
    datatype, idcol, firstcolpos,
    testfraction::Union{Float64,Nothing}=nothing, randseed=1892
)
    # read the file into a df
    if datatype == "Reads"
        df = DataFrame(CSV.File(infile, delim=delim))
    elseif datatype == "Proteins"
        df1 = DataFrame(CSV.File(infile, delim=delim, select=collect(1:firstcolpos-1)))
        # make sure columns of AAs containing only Ts aren't parsed as boolean columns
        df2 = DataFrame(CSV.File(infile, delim=delim, drop=collect(1:firstcolpos-1), types=String))
        df = hcat(df1, df2)
    else
        throw("datatype must be either `Reads` or `Proteins`")
    end

    # keep only rows edited at least once
    removeuneditedrows && !allrowsedited && filter!(row -> row[editingcol] > 0, df)
    # take a subset of the df for testing purposes
    if testfraction !== nothing
        Random.seed!(randseed)
        nrows = size(df, 1)
        nsamplerows = Int(round(testfraction * nrows))
        samplerows = sample(1:nrows, nsamplerows, replace=false)
        df = df[samplerows, :]
    end
    # flatten the df by exploding the `Reads` col, 
    # which denotes the reads supporting the unique observation (read / unique) the row represents
    transform!(df, :Reads => (x -> split.(x, ",")) => :Reads)
    df = flatten(df, :Reads)
    if datatype == "Reads"
        df = hcat(select(df, idcol), df[:, firstcolpos:end])
        firstcolpos = 2
    elseif datatype == "Proteins"
        df = hcat(select(df, idcol), toAAset.(df[:, firstcolpos:end]))
        firstcolpos = 2
    else
        throw("datatype must be either `Reads` or `Proteins`")
    end
    # remove uniformative cols (altough we apply the function on the idcol it shouldn't matter, as it is unique)
    df[:, map(col -> length(unique(col)) > 1, eachcol(df))]
    return df, firstcolpos
end


"""
The final result: 10 repetions of `f1` and `f2` over the neighborhood matrix of the graph, for each fraction of reads.  
Each such result contains a list of compatible unique reads (they are all different from each other).
"""
function main()
    # read & parse command-line args
    parsedargs = parsecmd()
    infiles = parsedargs["infiles"]
    delim = parsedargs["delim"]
    # samplesnames = parsedargs["samplesnames"]
    postfix = parsedargs["postfix"]
    idcol = parsedargs["idcol"]
    editingcol = parsedargs["editingcol"]
    firstcolpos = parsedargs["firstcolpos"]
    removeuneditedrows = parsedargs["removeuneditedrows"]
    allrowsedited = parsedargs["allrowsedited"]
    datatype = parsedargs["datatype"]
    outdir = parsedargs["outdir"]
    fracstep = parsedargs["fracstep"]
    maxfrac = parsedargs["maxfrac"]
    fracrepetitions = parsedargs["fracrepetitions"]
    algrepetitions = parsedargs["algrepetitions"]
    flatparallelism = parsedargs["flatparallelism"]
    testfraction = parsedargs["testfraction"]
    randseed = parsedargs["randseed"]

    # do the thing for each file, in ascending order of size
    samplesnames = map(x -> replace(splitdir(x)[2], postfix => ""), infiles)
    files_n_samples = [(infile, samplename) for (infile, samplename) ∈ zip(infiles, samplesnames)]
    sort!(files_n_samples, by=(x -> filesize(x[1])))
    infiles = map(x -> x[1], files_n_samples)
    samplesnames = map(x -> x[2], files_n_samples)

    # infiles = ["D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv"]
    # # protsfile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
    # delim = "\t"
    # samplesnames = ["GRIA-CNS-RESUB.C0x1291.MinRQ998"]
    # idcol = "Transcript"
    # # idcol = "Protein"
    # editingcol = "EditedPositions"
    # # editingcol = "MinNonSyns"
    # firstcolpos = 9
    # # firstcolpos = 15
    # removeuneditedrows = true
    # allrowsedited = false
    # datatype = "Reads"
    # # datatype = "Proteins"
    # outdir = "D.pealeii/MpileupAndTranscripts/RQ998.2"
    # fracstep = 1.0
    # maxfrac = 1.0
    # fracrepetitions = 1
    # algrepetitions = 1
    # flatparallelism = false

    for x ∈ eachindex(infiles)  # could be eachindex(samplesnames) as well
        infile = infiles[x]
        samplename = samplesnames[x]
        if flatparallelism
            flat_parallelism_submain(
                infile, delim, samplename, idcol, editingcol, firstcolpos, removeuneditedrows, allrowsedited, datatype,
                outdir, fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed
            )
        else
            nested_parallelism_submain(
                infile, delim, samplename, idcol, editingcol, firstcolpos, removeuneditedrows, allrowsedited, datatype,
                outdir, fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed
            )
        end
    end
end


function flat_parallelism_submain(
    infile, delim, samplename, idcol, editingcol, firstcolpos, removeuneditedrows, allrowsedited, datatype,
    outdir, fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed
)

    df, firstcolpos = preparedf(
        infile, delim,
        removeuneditedrows, allrowsedited, editingcol,
        datatype, idcol, firstcolpos,
        testfraction, randseed
    )

    # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
    G = indistinguishable_rows(df, idcol; firstcolpos)

    # having built G, we only need to keep the reads and unique reads/proteins they support
    select!(df, idcol)

    # define [fracstep, 2 fracstep, 3 fracstep, ..., maxfrac] fractions of the rows
    # (always include maxfrac in fractions, e.g., when fracstep=0.3 and maxfrac=1.0)
    fractions = collect(fracstep:fracstep:maxfrac) ∪ maxfrac
    # df of the results (each row represents a list of distinct unique read/protein under some conditions)
    results = emptyresults()

    nrows = size(df, 1)
    for fraction ∈ fractions

        nsamplerows = convert(Int, round(fraction * nrows))

        for fracrepetition ∈ 1:fracrepetitions  # repetitions for each fraction
            # assemble sub neighborhood lists of indistinguishable sampled rows by using the pre-computed complete graph
            sampleG = get_data_sample(G, fraction, nsamplerows, df, idcol)
            # find subsets of distinct unique samples (unique reads/proteins),using ascending & descending algorithms 
            # in $repetions in parallel
            distinctasc = [@spawn distinctascending(sampleG) for _ ∈ 1:algrepetitions]
            distinctdesc = [@spawn distinctdescending(sampleG) for _ ∈ 1:algrepetitions]
            distinctasc = [fetch(t) for t ∈ distinctasc] # unique samples obtained by the ascending algorithm
            distinctdesc = [fetch(t) for t ∈ distinctdesc] # unique samples obtained by the descending algorithm
            # add the compatible unique samples to the results
            results = emptyresults()
            for (algrepetition, distinctsamples) ∈ enumerate(distinctasc)
                push!(results, [fraction, fracrepetition, "Ascending", algrepetition, length(distinctsamples), distinctsamples])
            end
            for (algrepetition, distinctsamples) ∈ enumerate(distinctdesc)
                push!(results, [fraction, fracrepetition, "Descending", algrepetition, length(distinctsamples), distinctsamples])
            end

            # concatanate the vector of distinct unique samples in each result into a comma-separated string
            transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)
        end
    end
    # write the results to a csv file
    outfile = abspath(outdir) * "/" * samplename * ".DistinctUnique$datatype.csv"
    CSV.write(outfile, results; delim)
end


function nested_parallelism_submain(
    infile, delim, samplename, idcol, editingcol, firstcolpos, removeuneditedrows, allrowsedited, datatype,
    outdir, fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed
)

    df, firstcolpos = preparedf(
        infile, delim,
        removeuneditedrows, allrowsedited, editingcol,
        datatype, idcol, firstcolpos,
        testfraction, randseed
    )

    # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
    G = indistinguishable_rows(df, idcol; firstcolpos)

    # having built G, we only need to keep the reads and unique reads/proteins they support
    select!(df, idcol)

    # define [fracstep, 2 fracstep, 3 fracstep, ..., maxfrac] fractions of the rows
    # (always include maxfrac in fractions, e.g., when fracstep=0.3 and maxfrac=1.0)
    fractions = collect(fracstep:fracstep:maxfrac) ∪ maxfrac

    fractions_results = [
        @spawn run_fraction(df, idcol, G, fraction, fracrepetitions, algrepetitions)
        for fraction ∈ fractions
    ]
    fractions_results = [fetch(t) for t ∈ fractions_results]
    results = vcat(fractions_results...)

    # write the results to a csv file
    outfile = abspath(outdir) * "/" * samplename * ".DistinctUnique$datatype.csv"
    CSV.write(outfile, results; delim)
end


function run_fraction(
    df::DataFrame,
    idcol::String,
    G::Dict,
    fraction::Float64,
    fracrepetitions::Int64,
    algrepetitions::Int64
)
    nrows = size(df, 1)
    nsamplerows = convert(Int, round(fraction * nrows)) # number of rows to sample for this fraction of data    

    fracrepetitions_results = [
        @spawn run_fracrepetition(df, idcol, G, fraction, fracrepetition, algrepetitions, nsamplerows)
        for fracrepetition ∈ 1:fracrepetitions
    ]
    fracrepetitions_results = [fetch(t) for t ∈ fracrepetitions_results]

    results = vcat(fracrepetitions_results...)

    return results

end


function run_fracrepetition(
    df::DataFrame,
    idcol::String,
    G::Dict,
    fraction::Float64,
    fracrepetition::Int64,
    algrepetitions::Int64,
    nsamplerows::Int64
)
    # assemble sub neighborhood lists of indistinguishable sampled rows by using the pre-computed complete graph
    sampleG = get_data_sample(G, fraction, nsamplerows, df, idcol)

    # obtain $algrepetitions × 2 (the asc- and desc algorithms) results of distinct rows
    results = solve(sampleG, fraction, fracrepetition, algrepetitions)

    return results
end


emptyresults() = DataFrame(
    Fraction=Vector{Float64}(undef, 0),
    FractionRepetition=Vector{Int64}(undef, 0),
    Algorithm=Vector{String}(undef, 0),
    AlgorithmRepetition=Vector{Int64}(undef, 0),
    NumUniqueSamples=Vector{Int64}(undef, 0),
    UniqueSamples=Vector{Vector{String}}(undef, 0),
)


function get_data_sample(G::Dict, fraction::Float64, nsamplerows::Int64, df::DataFrame, idcol::String)
    if fraction < 1.0
        # sample a fraction of the rows (non-unique samples of reads/proteins)
        samplerows = sample(collect(1:size(df, 1)), nsamplerows, replace=false)
        # get their corresponding unique ids
        sampleids = unique(df[samplerows, idcol])
        # assemble sub neighborhood lists of uncompatible unique sampled rows by using the pre-computed complete graph
        sampleG = subgraph(G, sampleids)
    else  # fraction == 1.0, so there's no need to create a sub-graph
        sampleG = G
    end
    return sampleG
end


function solve(
    sampleG::Dict,
    fraction::Float64,
    fracrepetition::Int64,
    algrepetitions::Int64
)
    # find subsets of distinct unique samples (unique reads/proteins),using ascending & descending algorithms 
    # in $repetions in parallel
    distinctasc = [@spawn distinctascending(sampleG) for _ ∈ 1:algrepetitions]
    distinctdesc = [@spawn distinctdescending(sampleG) for _ ∈ 1:algrepetitions]
    distinctasc = [fetch(t) for t ∈ distinctasc] # unique samples obtained by the ascending algorithm
    distinctdesc = [fetch(t) for t ∈ distinctdesc] # unique samples obtained by the descending algorithm
    # add the compatible unique samples to the results
    results = emptyresults()
    for (algrepetition, distinctsamples) ∈ enumerate(distinctasc)
        push!(results, [fraction, fracrepetition, "Ascending", algrepetition, length(distinctsamples), distinctsamples])
    end
    for (algrepetition, distinctsamples) ∈ enumerate(distinctdesc)
        push!(results, [fraction, fracrepetition, "Descending", algrepetition, length(distinctsamples), distinctsamples])
    end
    # concatanate the vector of distinct unique samples in each result into a semicolon-separated string
    transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)
    # return results
    return results
end


main()



# @benchmark main()


# infile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv"
# samplename = "GRIA-CNS-RESUB.C0x1291.MinRQ998"
# firstcolpos = 9
# editingcol = "EditedPositions"
# idcol = "Transcript"
# datatype = "Reads"


# infile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
# samplename = "GRIA-CNS-RESUB.C0x1291.MinRQ998"
# firstcolpos = 15
# editingcol = "MinNonSyns"
# idcol = "Protein"
# allrowsedited = false
# removeuneditedrows = false
# outdir = "D.pealeii/MpileupAndTranscripts/RQ998.2"
# fracstep = 1.0
# maxfrac = 1.0
# fracrepetitions = 1
# algrepetitions = 1
# delim = "\t"
# datatype = "Proteins"
# fraction = 0.5
# fracrepetition = 1
# algrepetitions = 1



# df, firstcolpos = preparedf(
#     infile, delim,
#     removeunedited, allrowsedited, editingcol,
#     datatype, idcol, firstcolpos,
# )
# # read the file into a df
# df = DataFrame(CSV.File(infile, delim=delim))

# # keep only rows edited at least once
# if removeunedited && !allrowsedited
#     filter!(row -> row[editingcol] > 0, df)
# end

# transform!(df, :Reads => (x -> split.(x, ",")) => :Reads)
# df = flatten(df, :Reads)
# if datatype == "Reads"
#     df = hcat(select(df, idcol), df[:, firstcolpos:end])
#     firstcolpos = 2
# elseif datatype == "Proteins"
#     df = hcat(select(df, idcol), toAAset.(df[:, firstcolpos:end]))
#     firstcolpos = 2
# else
#     throw("datatype must be either `Reads` or `Proteins`")
# end

# df = df[vcat(collect(1:10), collect(size(df, 1)-100:size(df, 1))), :] # taking a subset for test purposes

# @time G = uncompatiblerows(df, idcol; firstcolpos)

# # having built G, we only need to keep the reads and unique reads/proteins they support
# select!(df, idcol)


# nrows = size(df, 1)
# nsamplerows = convert(Int, round(fraction * nrows)) # numbe

# # assemble sub neighborhood lists of uncompatible sampled rows by using the pre-computed complete graph
# sampleG = get_data_sample(G, fraction, nsamplerows, df, idcol)

# # # obtain $algrepetitions × 2 (the asc- and desc algorithms) results of compatible rows
# results = solve(sampleG, df, idcol, fraction, fracrepetition, algrepetitions)


# # find which rows are compatible with each other using ascending & descending algorithms in $repetions in parallel
# ascrows = [@spawn comptascending(sampleG) for _ ∈ 1:algrepetitions]
# descrows = [@spawn comptdescending(sampleG) for _ ∈ 1:algrepetitions]
# ascrows = [fetch(t) for t ∈ ascrows]
# descrows = [fetch(t) for t ∈ descrows]
# # extract the compatible unique samples (unique reads/proteins) by their rows' indices
# # uniqasc = [df[rows, idcol] for rows ∈ ascrows]  # unique samples obtained by the ascending algorithm
# # uniqdesc = [df[rows, idcol] for rows ∈ descrows]  # unique samples obtained by the descending algorithm
# uniqasc = ascrows  # unique samples obtained by the ascending algorithm
# uniqdesc = descrows  # unique samples obtained by the descending algorithm
# # add the compatible unique samples to the results
# results = emptyresults()
# for (algrepetition, uniqsamples) ∈ enumerate(uniqasc)
#     push!(results, [fraction, fracrepetition, "Ascending", algrepetition, length(uniqsamples), uniqsamples])
# end
# for (algrepetition, uniqsamples) ∈ enumerate(uniqdesc)
#     push!(results, [fraction, fracrepetition, "Descending", algrepetition, length(uniqsamples), uniqsamples])
# end
# # concatanate the vector of unique samples in each result into a semicolon-separated string
# transform!(results, :UniqueSamples => (x -> join.(x, ",")) => :UniqueSamples)


