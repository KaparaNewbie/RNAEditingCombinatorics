using DataFrames
using CSV
using Random # for MersenneTwister
using StatsBase  # for StatsBase.sample



function find_B_and_E_for_M(M)
    rows, cols = size(M)

    B = zeros(Int, rows)

    # E = Array{Int}(undef, rows)
    # E .= rows + 1
    E = fill(rows + 1, rows)

    # iterate over each row
    for i ∈ 1:rows
        # iterate from start to end
        for j ∈ 1:cols
            if M[i, j]
                B[i] = j
            else
                break
            end
        end
        # iterate from end to start
        for j ∈ reverse(1:cols)
            if M[i, j]
                E[i] = j
            else
                break
            end
        end
    end

    return B, E
end



function are_there_chimeric_reads(M)

    B, E = find_B_and_E_for_M(M)

    return maximum(B) >= minimum(E)

end


function find_all_chimeric_reads(M)

    B, E = find_B_and_E_for_M(M)
    
    sorted_B, idxs_B = sort(B), sortperm(B)
    sorted_E, idxs_E = sort(E), sortperm(E)

    i = 1 # for i in eachindex(sorted_B)
    k = searchsortedlast(sorted_E, sorted_B[i])

    pairs = Vector{Tuple{Int,Int}}()
    j = 1

    for i ∈ 1:length(sorted_B)
        # println("i: $i, j: $j")
        while j <= length(sorted_E) && sorted_E[j] <= sorted_B[i]
            j += 1
        end
        # println("j: $j")
        # All E[1:j-1] ≤ B[i]
        for k ∈ 1:j-1
            # println("k: $k")
            push!(pairs, (idxs_B[i], idxs_E[k]))
        end
    end

    return pairs

end



B = [7, 7, 3, 3, 0]
E = [4, 10, 4, 2, 11]











M = rand(Bool, (1000, 100))
# find_chimeric_reads(M)
are_there_chimeric_reads(M)


M = falses(1000, 100)
M[2, 1:45] .= true
M[3, 42:end] .= true
# @assert find_chimeric_reads(M2) == (2, 3)
# @assert find_chimeric_reads(M)
@assert are_there_chimeric_reads(M)


M = falses(1000, 100)
M[2, 1:45] .= true
M[3, 42:end] .= true
M[4, 40:end] .= true
@assert are_there_chimeric_reads(M)
find_all_chimeric_reads(M)

M = falses(1000, 100)
# @assert find_chimeric_reads(M) == nothing
# @assert !find_chimeric_reads(M)
@assert !are_there_chimeric_reads(M)
find_all_chimeric_reads(M)


function prepare_expression_df(expression_file)
    expression_df = DataFrame(CSV.File(expression_file, delim="\t"))

    # count the number of distinct proteins per solution
    expression_df = transform(groupby(expression_df, "#Solution"), nrow => "DistinctProteins")

    # filter the DataFrame to keep only the rows with the maximum number of distinct proteins
    expression_df = subset(expression_df, :DistinctProteins => x -> x .== maximum(x))

    # if there are multiple solutions with the same maximum number of distinct proteins, keep only one
    max_rand_solution = sample(MersenneTwister(1892), unique(expression_df[!, "#Solution"]), 1, replace=false)
    expression_df = subset(expression_df, "#Solution" => x -> x .== max_rand_solution)

    # sort the DataFrame by the total expression level, from the lowest to the highest
    expression_df = sort(expression_df, "TotalWeightedSupportingReads")

    # keep only the 1000 rarest and 300 most common distinct proteins
    # (the first 1000 rows and the last 300 rows)
    rare_expression_df = expression_df[1:1000, :]
    common_expression_df = expression_df[end-299:end, :]
    insertcols!(rare_expression_df, "ExpressionStatus" => "Rare")
    insertcols!(common_expression_df, "ExpressionStatus" => "Common")
    expression_df = vcat(rare_expression_df, common_expression_df)


    # expression_df[!, "Reads"] = replace.(expression_df[!, "Reads"], "SubString{String}" => "", "[\"" => "", "\"]" => "")
    # expression_df[!, "Reads"] = split.(expression_df[!, "Reads"], ",")
    # unique(length.(expression_df[!, "Reads"]))
    # expression_df[!, "AdditionalSupportingReadsIDs"]

    # filter out the columns we don't need
    cols_to_keep = split("""Gene
    Protein
    AdditionalSupportingProteinsIDs
    TotalWeightedSupportingReads
    ExpressionStatus""", "\n")
    expression_df = expression_df[:, cols_to_keep]

    expression_df[!, "AdditionalSupportingProteinsIDs"] = split.(expression_df[!, "AdditionalSupportingProteinsIDs"], ",")

    proteins = expression_df[:, "Protein"]
    additional_supporting_proteins = expression_df[:, "AdditionalSupportingProteinsIDs"]

    push!.(additional_supporting_proteins, proteins)

    insertcols!(expression_df, "OriginalAndAdditionalSupportingProteinsIDs" => additional_supporting_proteins)
    select!(expression_df, Not("AdditionalSupportingProteinsIDs"))

    expression_df = flatten(expression_df, "OriginalAndAdditionalSupportingProteinsIDs")
    transform!(expression_df, "OriginalAndAdditionalSupportingProteinsIDs" => "OriginalOrAdditionalSupportingProtein")
    select!(expression_df, Not("OriginalAndAdditionalSupportingProteinsIDs"))

    return expression_df
end


function prepare_unique_proteins_df(unique_proteins_file, firstcolpos)

    unique_proteins_df = DataFrame(CSV.File(unique_proteins_file, delim="\t", select=collect(1:firstcolpos-1), types=Dict("Protein" => String, "Reads" => String)))
    select!(unique_proteins_df, ["Gene", "Protein", "Transcripts"])
    unique_proteins_df[!, "Protein"] = InlineString.(unique_proteins_df[!, :Protein])
    transform!(unique_proteins_df, :Transcripts => (x -> split.(x, ",")) => :Transcripts)
    unique_proteins_df = flatten(unique_proteins_df, "Transcripts")
    transform!(unique_proteins_df, "Transcripts" => "UniqueRead")
    select!(unique_proteins_df, Not("Transcripts"))
    return unique_proteins_df

end


function prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)
    unique_reads_df = DataFrame(CSV.File(unique_reads_file, delim="\t", types=Dict("Reads" => String)))
    select!(unique_reads_df, vcat(["Gene", "Transcript"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
    rename!(unique_reads_df, "Transcript" => "UniqueRead")
    return unique_reads_df
end


expression_file = "TempData/GRIA.DistinctUniqueProteins.ExpressionLevels.csv"
unique_reads_file = "TempData/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz"
unique_proteins_file = "TempData/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
unique_proteins_first_col_pos = 15
unique_reads_first_col_pos = 9


expression_df = prepare_expression_df(expression_file)
unique_proteins_df = prepare_unique_proteins_df(unique_proteins_file, unique_proteins_first_col_pos)
unique_reads_df = prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)

unique_reads_and_proteins_df = innerjoin(unique_proteins_df, unique_reads_df, on=["Gene", "UniqueRead"])

expression_df = leftjoin(expression_df, unique_reads_and_proteins_df, on=["Gene", "OriginalOrAdditionalSupportingProtein" => "Protein"])

# now we have a dataframe with the original distinct proteins (300 most common and 1000 rarest) and their expression levels,
# as well as the indistinguishable proteins supporting them,
# and the unique reads supporting them
# (each row is a single unique read supporting a unique protein, wether it is a distinct or indistinguishable one which underwent reasassignment)



common_expression_df = expression_df[expression_df.ExpressionStatus .== "Common", :]
rare_expression_df = expression_df[expression_df.ExpressionStatus .== "Rare", :]

# original_common_proteins = unique(expression_df[expression_df.ExpressionStatus .== "Common", :Protein])
original_common_proteins = unique(common_expression_df[!, :Protein])

original_common_protein = original_common_proteins[1]
common_expression_df[common_expression_df.Protein .== original_common_protein, :]


unique(rare_expression_df)

