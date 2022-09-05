# using BioSequences
# using Transducers
# using LinearAlgebra



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
    @info "$(loggingtime())\tindistinguishable_rows"

    nrows = size(M, 1)
    length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )   # todo verify the integrity of parallelization with ThreadsX
    nodeseltype = eltype(ids)
    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row]) # todo verify the integrity of parallelization with ThreadsX
    emptyAAintersections = ThreadsX.Set(
        [(x, y)
         for x ∈ AAsets for y ∈ AAsets
         if x ∩ y == ∅]
    ) # todo verify the integrity of parallelization with ThreadsX
    Gs = tcollect(indistinguishable_vs_for_u(M, emptyAAintersections, ids, nodeseltype, i) for i ∈ 1:nrows)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
    return G
end


function indistinguishable_vs_for_u(
    M::Matrix{Set{AminoAcid}}, emptyAAintersections, ids, nodeseltype, i
)
    u = ids[i]
    @views rowᵢ = M[i, :]
    G = Dict(u => Set{nodeseltype}())
    nrows = size(M, 1)
    for j ∈ 1:nrows   # compare rowᵢ & rowⱼ if i != j   # todo use @incounds? https://docs.julialang.org/en/v1/base/base/#Base.@inbounds
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
    @info "$(loggingtime())\tindistinguishable_rows"

    # informativecols = tcollect(i for (i, col) in enumerate(eachcol(M)) if length(unique(col)) > 1)
    # M = M[:, informativecols]
    nrows = size(M, 1)
    length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    ) # todo verify the integrity of parallelization with ThreadsX
    nodeseltype = eltype(ids)
    @views M = replace(M, -1 => missing)
    Gs = tcollect(indistinguishable_vs_for_u(M, ids, nodeseltype, i) for i ∈ 1:nrows)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
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