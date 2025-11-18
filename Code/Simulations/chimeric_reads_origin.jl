using DataFrames
using CSV
using Random # for MersenneTwister
using StatsBase  # for StatsBase.sample
using CairoMakie
using Format


function find_B_and_E_for_M(M)
	rows, cols = size(M)

	B = zeros(Int, rows)

	# E = Array{Int}(undef, rows)
	# E .= rows + 1
	# E = fill(rows + 1, rows)
	E = fill(cols + 1, rows)

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


# r1 = M[1, :]
# r2 = M[2, :]

# r1_b = findfirst(x -> !x, r1) - 1
# ri_e = findlast(x -> !x, r1) + 1
# r2_b = findfirst(x -> !x, r2) - 1
# r2_e = findlast(x -> !x, r2) + 1

function find_all_chimeric_reads(M)

	B, E = find_B_and_E_for_M(M)

	sorted_B, idxs_B = sort(B), sortperm(B)
	sorted_E, idxs_E = sort(E), sortperm(E)

	i = 1 # for i in eachindex(sorted_B)
	k = searchsortedlast(sorted_E, sorted_B[i])

	pairs = Vector{Tuple{Int, Int}}()
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



# M = rand(Bool, (1000, 100))
# # find_chimeric_reads(M)
# are_there_chimeric_reads(M)


M = falses(1000, 100)
M[2, 1:45] .= true
M[3, 42:end] .= true
@assert are_there_chimeric_reads(M)
find_all_chimeric_reads(M)

M = falses(1000, 100)
M[2, 1:45] .= true
M[3, 42:end] .= true
M[4, 40:end] .= true
@assert are_there_chimeric_reads(M)
find_all_chimeric_reads(M)

M = falses(1000, 100)
@assert !are_there_chimeric_reads(M)
find_all_chimeric_reads(M)


function validate_chimeric_reads(rare_read_editing_status, two_common_reads_editing_statuses)

	# Validate the input dimensions:
	# rare_read_editing_status must be a vector of length m_sites,
	# two_common_reads_editing_statuses must have exactly two rows and m_sites columns.
	n_reads, m_sites = size(two_common_reads_editing_statuses)
	n_reads == 2 || error("two_common_reads_editing_statuses must have exactly two rows")
	if size(rare_read_editing_status) !== (m_sites,)
		rare_read_editing_status = rare_read_editing_status'
		size(rare_read_editing_status) !== (m_sites,) && error("rare_read_editing_status must be a vector of length $n_reads")
	end

	s1 = two_common_reads_editing_statuses[1, :] # first common read editing status
	s2 = two_common_reads_editing_statuses[2, :] # second common read editing status

	options = []

	for i ∈ eachindex(rare_read_editing_status) # could be one of the common reads as well

		# i = 101

		option1 = vcat(s1[begin:i-1], s2[i:end])
		option2 = vcat(s2[begin:i-1], s1[i:end])
		push!(options, option1)
		push!(options, option2)
	end

	return rare_read_editing_status in options

end



s1 = [1, 0, -1]
s2 = [-1, 0, 0]

s3 = [1, 0, 0]
s4 = [1, -1, 1]


@assert validate_chimeric_reads(s3, [s1 s2]')
@assert !(validate_chimeric_reads(s4, [s1 s2]'))


function prepare_expression_df(expression_file, x_common_proteins = 300, y_rare_proteins = 1000, discard_reassigned_reads::Bool = true)

	# TODO use an input largest solutuon rather than randomly selecting one

	expression_df = DataFrame(CSV.File(expression_file, delim = "\t"))

	# count the number of distinct proteins per solution
	expression_df = transform(groupby(expression_df, "#Solution"), nrow => "DistinctProteins")

	# filter the DataFrame to keep only the rows with the maximum number of distinct proteins
	expression_df = subset(expression_df, :DistinctProteins => x -> x .== maximum(x))

	# if there are multiple solutions with the same maximum number of distinct proteins, keep only one
	max_rand_solution = sample(MersenneTwister(1892), unique(expression_df[!, "#Solution"]), 1, replace = false)
	expression_df = subset(expression_df, "#Solution" => x -> x .== max_rand_solution)

	# sort the DataFrame by the total expression level, from the lowest to the highest
	expression_df = sort(expression_df, "TotalWeightedSupportingReads")

	# keep only the y_rare_proteins rarest and x_common_proteins most common distinct proteins
	# (the first y_rare_proteins rows and the last x_common_proteins rows)
	rare_expression_df = expression_df[1:y_rare_proteins, :]
	# common_expression_df = expression_df[end-x_common_proteins:end, :]
	common_expression_df = expression_df[end-x_common_proteins+1:end, :]
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

	# replace missing values in "AdditionalSupportingProteinsIDs" with empty strings
	expression_df[!, "AdditionalSupportingProteinsIDs"] = coalesce.(expression_df[!, "AdditionalSupportingProteinsIDs"], "")

	# split the "AdditionalSupportingProteinsIDs" column into arrays
	expression_df[!, "AdditionalSupportingProteinsIDs"] = split.(expression_df[!, "AdditionalSupportingProteinsIDs"], ",")

	proteins = expression_df[:, "Protein"]
	additional_supporting_proteins = expression_df[:, "AdditionalSupportingProteinsIDs"]

	push!.(additional_supporting_proteins, proteins)
	# when pushing the protein to the additional_supporting_proteins, 
	# after replacing some missing values with empty strings,
	# we have some arrays with empty values (which are the result of pushing the original protein to the 
	# missing supporting proteins), e.g.:
	# 1300-element Vector{Vector{SubString{String}}}:
	# ["", "vz"]
	# ["", "wr"]
	additional_supporting_proteins = filter.(!isempty, additional_supporting_proteins)

	insertcols!(expression_df, "OriginalAndAdditionalSupportingProteinsIDs" => additional_supporting_proteins)
	select!(expression_df, Not("AdditionalSupportingProteinsIDs"))

	expression_df = flatten(expression_df, "OriginalAndAdditionalSupportingProteinsIDs")
	transform!(expression_df, "OriginalAndAdditionalSupportingProteinsIDs" => "OriginalOrAdditionalSupportingProtein")
	select!(expression_df, Not("OriginalAndAdditionalSupportingProteinsIDs"))

	# if we want to discard the reassigned reads, we filter out the rows where "Protein" != "OriginalOrAdditionalSupportingProtein"
	if discard_reassigned_reads
		expression_df = expression_df[expression_df.Protein.==expression_df.OriginalOrAdditionalSupportingProtein, :]
	end

	return expression_df
end


function softcomparison(r1, r2)
	(size(r1) == size(r2) && length(size(r1)) == length(size(r2)) == 1) || error("r1 and r2 must have the same size")

	# initialize a boolean array of the same size as r1 and r2 with all true values
	result = trues(size(r1))

	# only a strong disagreement (0 vs 1) makes a difference, 
	# while a missing value (-1) is compatible with both 0 and 1
	for i ∈ eachindex(result)
		if r1[i] != -1 && r2[i] != -1 && r1[i] != r2[i]
			result[i] = false
		end
	end

	return result
end


function prepare_unique_proteins_df(unique_proteins_file, firstcolpos)

	unique_proteins_df = DataFrame(CSV.File(unique_proteins_file, delim = "\t", select = collect(1:firstcolpos-1), types = Dict("Protein" => String, "Reads" => String)))
	select!(unique_proteins_df, ["Gene", "Protein", "Transcripts"])
	unique_proteins_df[!, "Protein"] = InlineString.(unique_proteins_df[!, :Protein])
	transform!(unique_proteins_df, :Transcripts => (x -> split.(x, ",")) => :Transcripts)
	unique_proteins_df = flatten(unique_proteins_df, "Transcripts")
	transform!(unique_proteins_df, "Transcripts" => "UniqueRead")
	select!(unique_proteins_df, Not("Transcripts"))
	return unique_proteins_df

end


function prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)
	unique_reads_df = DataFrame(CSV.File(unique_reads_file, delim = "\t", types = Dict("Reads" => String)))
	select!(unique_reads_df, vcat(["Gene", "Transcript"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	rename!(unique_reads_df, "Transcript" => "UniqueRead")
	return unique_reads_df
end




struct ReassignmentMetadata
	Protein::Any
	OriginalOrAdditionalSupportingProtein::Any
	TotalWeightedSupportingReads::Any
end


function process_one_sample(
	sample_name,
	expression_file,
	unique_proteins_file, unique_proteins_first_col_pos,
	unique_reads_file, unique_reads_first_col_pos,
	x_common_proteins = 300, y_rare_proteins = 1000,
	discard_reassigned_reads::Bool = true,
	soft_comparison::Bool = false,
)

	expression_df = prepare_expression_df(expression_file, x_common_proteins, y_rare_proteins, discard_reassigned_reads)
	unique_proteins_df = prepare_unique_proteins_df(unique_proteins_file, unique_proteins_first_col_pos)
	unique_reads_df = prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)

	unique_reads_and_proteins_df = innerjoin(unique_proteins_df, unique_reads_df, on = ["Gene", "UniqueRead"])

	expression_df = leftjoin(expression_df, unique_reads_and_proteins_df, on = ["Gene", "OriginalOrAdditionalSupportingProtein" => "Protein"])

	# now we have a dataframe with the original distinct proteins (300 most common and 1000 rarest) and their expression levels,
	# as well as the indistinguishable proteins supporting them,
	# and the unique reads supporting them
	# (each row is a single unique read supporting a unique protein, wether it is a distinct or indistinguishable one which underwent reasassignment)

	common_expression_df = expression_df[expression_df.ExpressionStatus.=="Common", :]
	rare_expression_df = expression_df[expression_df.ExpressionStatus.=="Rare", :]

	new_first_reads_col_pos = 7

	# common_expression_df[:, begin:new_first_reads_col_pos-1]

	gdf = groupby(common_expression_df[:, ["Gene", "ExpressionStatus", "UniqueRead", "Protein", "OriginalOrAdditionalSupportingProtein", "TotalWeightedSupportingReads"]], ["Gene", "ExpressionStatus", "UniqueRead"])

	# unique(common_expression_df[!, "UniqueRead"])
	# length(gdf)

	# struct ReassignmentMetadata
	#     Protein
	#     OriginalOrAdditionalSupportingProtein
	#     TotalWeightedSupportingReads
	# end


	compact_common_expression_df = combine(gdf) do subdf
		(; Metadata = hcat(collect(ReassignmentMetadata.(subdf.Protein, subdf.OriginalOrAdditionalSupportingProtein, subdf.TotalWeightedSupportingReads))))
	end

	# compact_common_expression_df[!, "Metadata"][1]
	# compact_common_expression_df[!, "Metadata"][1][1]

	# unique(common_expression_df[!, new_first_reads_col_pos-1:end], "UniqueRead")

	compact_common_expression_df = leftjoin(
		compact_common_expression_df,
		unique(common_expression_df[!, new_first_reads_col_pos-1:end], "UniqueRead"),
		on = ["UniqueRead"],
	)

	compact_common_expression_first_reads_col_pos = 5

	# common_unique_reads_editing_status = compact_common_expression_df[:, compact_common_expression_first_reads_col_pos:end]

	results = []

	original_rare_proteins = unique(rare_expression_df[!, :Protein])

	# TODO uncomment
	# one_original_rare_protein = original_rare_proteins[1]

	for one_original_rare_protein ∈ original_rare_proteins

		unique_reads_of_one_original_rare_protein_expression_df = rare_expression_df[rare_expression_df.Protein.==one_original_rare_protein, :]
		unique_reads_of_one_original_rare_protein_expression = unique(unique_reads_of_one_original_rare_protein_expression_df[!, "UniqueRead"])

		# TODO uncomment
		# one_unique_read_of_one_original_rare_protein = unique_reads_of_one_original_rare_protein_expression[1]

		for one_unique_read_of_one_original_rare_protein ∈ unique_reads_of_one_original_rare_protein_expression

			one_unique_read_of_one_original_rare_protein_df = unique_reads_of_one_original_rare_protein_expression_df[unique_reads_of_one_original_rare_protein_expression_df.UniqueRead.==one_unique_read_of_one_original_rare_protein, :][1, :] # take the first row of the unique reads of the original rare protein
			one_unique_read_of_one_original_rare_protein_editing_status = one_unique_read_of_one_original_rare_protein_df[new_first_reads_col_pos:end]


			common_unique_reads_editing_status = compact_common_expression_df[compact_common_expression_df.UniqueRead.!==one_unique_read_of_one_original_rare_protein, compact_common_expression_first_reads_col_pos:end]

			# one_unique_read_of_one_original_rare_protein_editing_status .== common_unique_reads_editing_status

			# common_unique_reads_editing_status .== one_unique_read_of_one_original_rare_protein_editing_status'

			one_unique_read_of_one_original_rare_protein_editing_status_array = Array(one_unique_read_of_one_original_rare_protein_editing_status)
			common_unique_reads_editing_status_array = Array(common_unique_reads_editing_status)

			# Compare the rare read to all common reads (row-wise)
			if soft_comparison
				# r1 = common_unique_reads_editing_status_array[1, end-9:end-2]
				# r2 = one_unique_read_of_one_original_rare_protein_editing_status_array[end-9:end-2]

				# r1 = collect(eachrow(common_unique_reads_editing_status_array))[1]
				# r2 = one_unique_read_of_one_original_rare_protein_editing_status_array

				M = softcomparison.(eachrow(common_unique_reads_editing_status_array), Ref(one_unique_read_of_one_original_rare_protein_editing_status_array))
				# M2 = softcomparison.(eachrow(common_unique_reads_editing_status_array[1:3, end-9:end-2]), Ref(one_unique_read_of_one_original_rare_protein_editing_status_array[end-9:end-2]))
				M = hcat(M...)'
			else
				M = common_unique_reads_editing_status_array .== one_unique_read_of_one_original_rare_protein_editing_status_array'
			end


			# size(M)[2]

			# sum(sum(row) == size(M)[2] for row ∈ eachrow(M)) # count the number of common reads that are identical to the rare read
			# sum(M, dims=2)

			is_chimeric = are_there_chimeric_reads(M)
			all_chimeric_reads_indices = find_all_chimeric_reads(M)

			result = (
				Sample = sample_name,
				XCommonProteins = x_common_proteins,
				YRareProteins = y_rare_proteins,
				ReassignendReadsDiscarded = discard_reassigned_reads,
				Protein = one_original_rare_protein,
				UniqueRead = one_unique_read_of_one_original_rare_protein,
				IsChimeric = is_chimeric,
				ChimericReadsIndices = all_chimeric_reads_indices,
				NumOfChimericCombinations = length(all_chimeric_reads_indices),
			)
			push!(results, result)

		end

	end

	results_df = DataFrame(results)

	return results_df

end





function per_platform_stats_df(
	X_common_proteins,
	Y_rare_proteins,
	platform,
	samples,
	unique_reads_files,
	unique_proteins_files,
	expression_files,
	unique_reads_first_col_pos,
	unique_proteins_first_col_pos,
	discard_reassigned_reads::Bool = true,
)
	x_and_y_common_and_rare_proteins = [(x, y) for x in X_common_proteins for y in Y_rare_proteins]

	results_dfs = []

	for soft_comparison in [true, false]

		for (x_common_proteins, y_rare_proteins) ∈ x_and_y_common_and_rare_proteins

			println("Processing samples with $x_common_proteins common proteins and $y_rare_proteins rare proteins...")

			for (sample_name, expression_file, unique_reads_file, unique_proteins_file) ∈ zip(
				samples, expression_files, unique_reads_files, unique_proteins_files)

				println("Processing sample: $sample_name")
				results_df = process_one_sample(
					sample_name,
					expression_file,
					unique_proteins_file, unique_proteins_first_col_pos,
					unique_reads_file, unique_reads_first_col_pos,
					x_common_proteins, y_rare_proteins,
					discard_reassigned_reads,
					soft_comparison,
				)
				insertcols!(results_df, "IsSoftComparison" => soft_comparison, "Platform" => platform)

				chimeric_results_df = results_df[results_df.IsChimeric.==true, :]
				println("Number of chimeric reads in $sample_name: $(nrow(chimeric_results_df))")
				# println(chimeric_results_df)
				push!(results_dfs, results_df)
			end
		end

	end


	results_df = vcat(results_dfs...)

	# grouped_results_df = groupby(results_df, [:Sample,])
	# combine(grouped_results_df, :IsChimeric => sum => :NumOfChimericReads,)

	# combine(groupby(results_df, [:Sample, :XCommonProteins, :YRareProteins, :ReassignendReadsDiscarded]), :IsChimeric => sum => :NumOfChimericReads)
	# combine(groupby(results_df, [:ReassignendReadsDiscarded, :XCommonProteins, :YRareProteins, :Sample,  ]), :IsChimeric => sum => :NumOfChimericReads)
	stats_df = combine(groupby(results_df, [:Platform, :IsSoftComparison, :XCommonProteins, :YRareProteins, :Sample]), :IsChimeric => sum => :NumOfChimericReads)
	insertcols!(stats_df, "%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.YRareProteins)

end






X_common_proteins = [300, 500, 1_000]
Y_rare_proteins = [1_000, 5_000, 10_000]



platform = "PacBio"
samples = ["GRIA2", "PCLO", "ADAR1", "IQEC1"]
unique_reads_files = [
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
]
unique_proteins_files = [
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
]
expression_files = [
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
]
unique_reads_first_col_pos = 9
unique_proteins_first_col_pos = 15

pacbio_stats_df = per_platform_stats_df(
	X_common_proteins,
	Y_rare_proteins,
	platform,
	samples,
	unique_reads_files,
	unique_proteins_files,
	expression_files,
	unique_reads_first_col_pos,
	unique_proteins_first_col_pos,
)







# TODO uncomment to test one sample
# sample_name = samples[1]
# expression_file = expression_files[1]
# unique_reads_file = unique_reads_files[1]
# unique_proteins_file = unique_proteins_files[1]
# x_common_proteins = 300
# y_rare_proteins = 1000

# results_df = process_one_sample(
#     sample_name,    
#     expression_file, 
#     unique_proteins_file, unique_proteins_first_col_pos,
#     unique_reads_file, unique_reads_first_col_pos
#     )

# results_df[results_df.IsChimeric .== true, :]


# describe(results_df[results_df.IsChimeric .== true, :NumOfChimericCombinations])

# x_and_y_common_and_rare_proteins = [
#     (300, 1_000), # x_common_proteins, y_rare_proteins
#     (1_000, 10_000),
# ]



function plot_stats_df(
	X_common_proteins,
	Y_rare_proteins,
	stats_df,
	samples,
	yscale = identity,
	colors = Makie.wong_colors()
)
	samples_cat = Dict(sample => i for (i, sample) in enumerate(samples))
	xticks = (
		1:length(samples),
		samples,
	)

	for soft_comparison in [true, false]

		fig = Figure(size = (600, 600))
		axes = []
		for (i, x_common_proteins) in enumerate(X_common_proteins)
			for (j, y_rare_proteins) in enumerate(Y_rare_proteins)
				title = "Common: $(format(x_common_proteins, commas=true))\nRare: $(format(y_rare_proteins, commas=true))"
				ax = Axis(
					fig[i, j],
					xticks = xticks,
					xticklabelrotation=π/4,
					subtitle = title,
					yscale = yscale
				)
				push!(axes, ax)
				subdf = stats_df[stats_df.XCommonProteins.==x_common_proteins.&&stats_df.YRareProteins.==y_rare_proteins.&&stats_df.IsSoftComparison.==soft_comparison, :]
				xs = [samples_cat[sample] for sample in subdf.Sample]
				ys = subdf[!, "%OfChimericReads"]

				barplot!(xs, ys,
					color = colors[xs],
					# strokecolor = :black, 
					# strokewidth = 1
				)
			end
		end
		linkxaxes!(axes...)
		linkyaxes!(axes...)
		Label(fig[begin:end, 0], "% of chimeric reads", rotation = pi / 2)  # y axis title

		if soft_comparison
			title = "Soft comparison"
		else
			title = "Strict comparison"
		end
		Label(fig[0, begin+1:end], title, fontsize = 18)  # main title

		display(fig)
	end

end




plot_stats_df(
	X_common_proteins,
	Y_rare_proteins,
	pacbio_stats_df,
	samples,
	log2
)


pacbio_stats_df[pacbio_stats_df.Sample .== "ADAR1" .&& pacbio_stats_df.IsSoftComparison, :]
pacbio_stats_df[pacbio_stats_df.Sample .== "ADAR1" .&& pacbio_stats_df.IsSoftComparison .== false, :]

pacbio_stats_df[pacbio_stats_df.Sample .== "IQEC1" .&& pacbio_stats_df.IsSoftComparison, :]
pacbio_stats_df[pacbio_stats_df.Sample .== "IQEC1" .&& pacbio_stats_df.IsSoftComparison .== false, :]






# (11, 17803)
# compact_common_expression_df[[11, 17803], 5:end]
# # sum([4020 in p for p in all_chimeric_reads_indices]) # 4020 is the index of the first common read

# # [p for p in all_chimeric_reads_indices if !(4020 in p)] # find all pairs with the first common read

# using DataStructures

# # counter([1, 2, 3])
# # unique.([[1, 2], [1, 1], [2, 3], [1, 2], [1, 3]])
# # vcat([[1, 2], [1]]...)
# # vcat(unique.([[1, 2], [1, 1], [2, 3], [1, 2], [1, 3]])...)
# # counter(vcat(unique.([[1, 2], [1, 1], [2, 3], [1, 2], [1, 3]])...))


# unique_chimeric_reads_indices_counter = counter(vcat(unique.(all_chimeric_reads_indices)...))
# # sort(unique_chimeric_reads_indices_counter)

# # unique_chimeric_reads_indices_counter[4020]

# unique_chimeric_reads_indices_counter_df = DataFrame(Indices=collect(keys(unique_chimeric_reads_indices_counter)), Counts=collect(values(unique_chimeric_reads_indices_counter)))
# sort!(unique_chimeric_reads_indices_counter_df, :Counts, rev=true)
# unique_chimeric_reads_indices_counter_df[unique_chimeric_reads_indices_counter_df.Counts .>= 1000, :]




# two_common_reads_editing_statuses = common_unique_reads_editing_status_array[[11, 17803], :]

# rare_read_editing_status = one_unique_read_of_one_original_rare_protein_editing_status_array'

# v = vcat(two_common_reads_editing_statuses, rare_read_editing_status)

# first_rare_read_editing_status



# validate_chimeric_reads(rare_read_editing_status, two_common_reads_editing_statuses)


# compact_common_expression_df[4020, :]
# compact_common_expression_df[4020, "Metadata"]


# sum(Array(compact_common_expression_df[4020, compact_common_expression_first_reads_col_pos:end]) .== one_unique_read_of_one_original_rare_protein_editing_status_array)

# compact_common_expression_df[17355, :] # this is the second common read that is chimeric with the rare read

# sum(Array(compact_common_expression_df[17355, compact_common_expression_first_reads_col_pos:end]) .== one_unique_read_of_one_original_rare_protein_editing_status_array)


# # compact_common_expression_df[[all_chimeric_reads_indices[1]...], :]
# two_common_reads_editing_statuses = common_unique_reads_editing_status_array[[all_chimeric_reads_indices[1]...], :]

# rare_read_editing_status = one_unique_read_of_one_original_rare_protein_editing_status_array'



# validate_chimeric_reads(rare_read_editing_status, two_common_reads_editing_statuses)

# M = two_common_reads_editing_statuses .== rare_read_editing_status

# are_there_chimeric_reads(M)

# find_all_chimeric_reads(M)