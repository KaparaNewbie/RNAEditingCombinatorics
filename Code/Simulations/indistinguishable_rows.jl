include(joinpath(@__DIR__, "issimilar.jl")) # for anysimilarity




"""
    indistinguishable_rows(M, ids)

Construct a simple graph `G` represented by a dictionary of neighborhood lists.  
Each row of `M` (a unique read) is a vertex in the graph.  
Two unique reads `i` and `j`, `i ≠ j`, are connected (`j ∈ G[i]`) if they cannot be distinguished by any editing position 
(the first of which is the column whose index is `firstcolpos`). 
That is, `i` and `j` are uncompatible and thus cannot be considered both as unique reads in final analysis. 
Else, `j ∉ G[i]`.
"""
function indistinguishable_rows(
    M::Matrix{Set{AminoAcid}}, ids,
    aagroups::Dict{AminoAcid,String}
)
    @info "$(loggingtime())\tindistinguishable_rows" aagroups

    nrows = size(M, 1)
    length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )
    nodeseltype = eltype(ids)
    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
    distinctAAsets = ThreadsX.Set(
        [
        (x, y)
        for x ∈ AAsets for y ∈ AAsets
        if !anysimilarity(x, y, aagroups)
    ]
    )
    Gs = tcollect(indistinguishable_vs_for_u(M, distinctAAsets, ids, nodeseltype, i) for i ∈ 1:nrows)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
    return G
end




"""
    indistinguishable_rows(M, ids)

Construct a simple graph `G` represented by a dictionary of neighborhood lists.  
Each row of `M` (a unique read) is a vertex in the graph.  
Two unique reads `i` and `j`, `i ≠ j`, are connected (`j ∈ G[i]`) if they cannot be distinguished by any editing position 
(the first of which is the column whose index is `firstcolpos`). 
That is, `i` and `j` are uncompatible and thus cannot be considered both as unique reads in final analysis. 
Else, `j ∉ G[i]`.
"""
function indistinguishable_rows(
    M::Matrix{Set{AminoAcid}}, ids,
    substitutionmatrix::SubstitutionMatrix{AminoAcid,Int64}, similarityscorecutoff::Int64, similarityvalidator::Function
)
    @info "$(loggingtime())\tindistinguishable_rows" substitutionmatrix similarityscorecutoff similarityvalidator

    nrows = size(M, 1)
    length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )
    nodeseltype = eltype(ids)
    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
    distinctAAsets = ThreadsX.Set(
        [
        (x, y)
        for x ∈ AAsets for y ∈ AAsets
        if !anysimilarity(x, y, substitutionmatrix, similarityscorecutoff, similarityvalidator)
    ]
    )
    Gs = tcollect(indistinguishable_vs_for_u(M, distinctAAsets, ids, nodeseltype, i) for i ∈ 1:nrows)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
    return G
end





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
    )
    nodeseltype = eltype(ids)

    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
    # sets whose intersection is empty
    distinctAAsets = ThreadsX.Set(
        [(x, y)
         for x ∈ AAsets for y ∈ AAsets
         if x ∩ y == ∅]
    )

    Gs = tcollect(indistinguishable_vs_for_u(M, distinctAAsets, ids, nodeseltype, i) for i ∈ 1:nrows)
    G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
    return G
end


function indistinguishable_vs_for_u(
    M::Matrix{Set{AminoAcid}},
    distinctAAsets,
    ids, nodeseltype, i
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
            # if (x, y) ∈ emptyAAintersections
            if (x, y) ∈ distinctAAsets
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

    nrows = size(M, 1)
    length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
        "Each id (a vertex in the graph) should be unique. " *
        "There should be a bijection between M's rows and ids, in order of appearence."
    )
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



function indistinguishable_rows(
    df::DataFrame, idcol;
    firstcolpos::Int=2, areuniuqe::Bool=false
)
    M, ids = preprocess(df, idcol, firstcolpos; areuniuqe)
    return indistinguishable_rows(M, ids)
end



function indistinguishable_rows(
    df::DataFrame, idcol,
    substitutionmatrix::SubstitutionMatrix{AminoAcid,Int64}, similarityscorecutoff::Int64, similarityvalidator::Function;
    firstcolpos::Int=2, areuniuqe::Bool=false
)
    M, ids = preprocess(df, idcol, firstcolpos; areuniuqe)
    return indistinguishable_rows(M, ids, substitutionmatrix, similarityscorecutoff, similarityvalidator)
end


function indistinguishable_rows(
    df::DataFrame, idcol,
    aagroups::Dict{AminoAcid,String};
    firstcolpos::Int=2, areuniuqe::Bool=false
)
    M, ids = preprocess(df, idcol, firstcolpos; areuniuqe)
    return indistinguishable_rows(M, ids, aagroups)
end




# todo find a more descriptive name for this function
function preprocess(df, idcol, firstcolpos; areuniuqe::Bool=false)
    if !areuniuqe
        df = unique(df, idcol)
        areuniuqe = true
    end
    M = Matrix(df[:, firstcolpos:end])
    ids = df[:, idcol]
    return M, ids
end