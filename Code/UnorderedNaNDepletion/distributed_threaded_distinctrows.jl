if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
else
    using BenchmarkTools
end


using BioSequences # for toAAset
using CSV # for CSV.read
using Random  # for Random.seed!(1892)
using StatsBase  # for StatsBase.sample
using DataFrames
using Transducers  # for tcollect


toAAset(x) = Set(map(aa -> convert(AminoAcid, only(aa)), split(x, ",")))


function preparedf(
    infile, delim,
    removeuneditedrows, allrowsedited, editingcol,
    datatype, idcol, firstcolpos,
    testfraction::Union{Float64,Nothing}=nothing, randseed=1892
)
    # read the file into a df
    df = DataFrame(CSV.File(infile, delim=delim))
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
    return df, firstcolpos
end


const ∅ = Set()


"""
    uncompatiblerows(M, firstcolpos)

Construct a simple graph `G` represented by a dictionary of neighborhood lists.  
Each row of `M` (a unique read) is a vertex in the graph.  
Two unique reads `i` and `j`, `i ≠ j`, are connected (`j ∈ G[i]`) if they cannot be distinguished by any of their features 
(the first of which is the column whose index is `firstcolpos`). 
That is, `i` and `j` are indistinguishable and thus cannot be considered both as unique samples in final analysis. 
Else, `j ∉ G[i]`.
"""


function sequential_indistinguishable_rows(M::Union{Matrix{Float64},Matrix{Int}}, ids; areuniuqe::Bool=false)
    if !areuniuqe
        @views M = unique(M, dims=1)
        ids = unique(ids)
    end
    nrows = size(M, 1)
    @assert length(ids) == nrows
    nodeseltype = eltype(ids)
    G = Dict(u => Set{nodeseltype}() for u ∈ ids)
    @views M = replace(M, -1 => missing)
    for i ∈ 1:nrows
        @views rowᵢ = M[i, :]
        u = ids[i]
        for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
            i == j && continue
            @views rowⱼ = M[j, :]
            v = ids[j]
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
    end
    return G
end

function sequential_indistinguishable_rows(M::Matrix{Set{AminoAcid}}, ids; areuniuqe::Bool=false)
    if !areuniuqe
        @views M = unique(M, dims=1)
        ids = unique(ids)
    end
    nrows = size(M, 1)
    @assert length(ids) == nrows
    nodeseltype = eltype(ids)
    G = Dict(u => Set{nodeseltype}() for u ∈ ids)
    AAsets = Set([x for row ∈ eachrow(M) for x ∈ row])
    emptyAAintersections = Set(
        [(x, y)
         for x ∈ AAsets for y ∈ AAsets
         if x ∩ y == ∅]
    )
    for i ∈ 1:nrows
        @views rowᵢ = M[i, :]
        u = ids[i]
        for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j
            i == j && continue
            @views rowⱼ = M[j, :]
            v = ids[j]
            distinct = false
            for (x, y) ∈ zip(rowᵢ, rowⱼ)
                # if we can find at least one col for which the rows' intersection is empty,
                # that is, if the rows don't share any possible amino acid for that position,
                # then we can consider these two rows as different (compatible)
                if (x, y) ∈ emptyAAintersections
                    distinct = true
                    break
                end
            end
            # if we can't differ the rows, we consider them as indistinguishable one from each other, 
            # and thus add them to the neighborhood graph
            !distinct && push!(G[u], v)
        end
    end
    return G
end


function threaded_indistinguishable_rows(M::Matrix{Set{AminoAcid}}, ids; areuniuqe::Bool=false)
    if !areuniuqe
        @views M = unique(M, dims=1)
        ids = unique(ids)
    end
    nrows = size(M, 1)
    @assert length(ids) == nrows
    nodeseltype = eltype(ids)
    AAsets = Set([x for row ∈ eachrow(M) for x ∈ row])
    emptyAAintersections = Set(
        [(x, y)
         for x ∈ AAsets for y ∈ AAsets
         if x ∩ y == ∅]
    )
    Gs = tcollect(indistinguishable_vs_for_u(M, emptyAAintersections, ids, nodeseltype, i) for i ∈ 1:nrows)
    # G = merge(Gs...)
    G = Dict(k => v for g in Gs for (k, v) in g) # each g has a single (k, v) pair
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


function threaded_indistinguishable_rows(M::Union{Matrix{Float64},Matrix{Int}}, ids; areuniuqe::Bool=false)
    if !areuniuqe
        @views M = unique(M, dims=1)
        ids = unique(ids)
    end
    nrows = size(M, 1)
    @assert length(ids) == nrows
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



function indistinguishable_rows(df::DataFrame, idcol, method::Symbol; firstcolpos::Int=2, areuniuqe::Bool=false)
    if !areuniuqe
        df = unique(df, idcol)
        areuniuqe = true
    end
    M = Matrix(df[:, firstcolpos:end])
    ids = df[:, idcol]
    if method == :sequential
        return sequential_indistinguishable_rows(M, ids; areuniuqe)
    elseif method == :threaded
        return threaded_indistinguishable_rows(M, ids; areuniuqe)
    else
        error("unknown method: " + method)
    end
end


function verify_G(G, df, idcol)
    df = unique(df, idcol)
    M = Matrix(df[:, firstcolpos:end])
    ids = df[:, idcol]
    # M = unique(M, dims=1)
    # ids = unique(ids)
    @views M = replace(M, -1 => missing)
    Gmissing = Dict()
    Gunwanted = Dict()
    nrows = size(M, 1)
    for (u, i) in zip(ids, 1:nrows)
        @views rowᵢ = M[i, :]
        present_vs = G[u]
        missing_vs = Set()
        unwanted_vs = Set()
        for (v, j) in zip(ids, 1:nrows)
            if i == j
                continue
            end
            @views rowⱼ = M[j, :]
            distinct = false
            for (x, y) ∈ zip(rowᵢ, rowⱼ)
                if x !== missing && y !== missing && x != y
                    distinct = true
                    break
                end
            end

            # !distinct && !(v in present_vs) && print("missing\n$u = $rowᵢ\n$v = $rowⱼ\n\n")
            # distinct && v in present_vs && print("unwanted\n$u = $rowᵢ\n$v = $rowⱼ\n\n")

            !distinct && !(v in present_vs) && push!(missing_vs, v)
            distinct && v in present_vs && push!(unwanted_vs, v)
        end
        Gmissing[u] = missing_vs
        Gunwanted[u] = unwanted_vs
    end
    return Gmissing, Gunwanted
end




# infile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv"
# idcol = "Transcript"
# editingcol = "EditedPositions"
# firstcolpos = 9
# datatype = "Reads"

infile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
idcol = "Protein"
editingcol = "MinNonSyns"
firstcolpos = 15
datatype = "Proteins"

delim = "\t"
removeuneditedrows = false
allrowsedited = false
# testfraction = 0.1
testfraction = nothing
randseed = 1892


df, firstcolpos = preparedf(
    infile, delim,
    removeuneditedrows, allrowsedited, editingcol,
    datatype, idcol, firstcolpos,
    testfraction, randseed
)



Gs = indistinguishable_rows(df, idcol, :sequential; firstcolpos)
Gs_missing, Gs_unwanted = verify_G(Gs, df, idcol);
for (k, v) in Gs_missing
    v != Set() && println("$k, $v")
end
for (k, v) in Gs_unwanted
    v != Set() && println("$k, $v")
end



Gt = indistinguishable_rows(df, idcol, :threaded; firstcolpos)
Gt_missing, Gt_unwanted = verify_G(Gt, df, idcol);
for (k, v) in Gt_missing
    v != Set() && println("$k, $v")
end
for (k, v) in Gt_unwanted
    v != Set() && println("$k, $v")
end

# proof that the threaded version is safe
tests = 100
Gts = [indistinguishable_rows(df, idcol, :threaded; firstcolpos) for _ in 1:tests]
for i in 1:tests-1
    @assert Gts[i] == Gts[i+1]
end

# proof that the threaded version gives the same result as the sequential one
Gs == Gt

# # benchmarks
@benchmark indistinguishable_rows(df, idcol, :threaded; firstcolpos) # threaded
@benchmark indistinguishable_rows(df, idcol, :sequential; firstcolpos) # sequential



