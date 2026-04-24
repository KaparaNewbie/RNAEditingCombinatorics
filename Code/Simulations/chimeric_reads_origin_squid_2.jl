using DataFrames
using CSV
using Random # for MersenneTwister
using StatsBase  # for StatsBase.sample
using CairoMakie
using Format


function find_B_and_E_for_M(M)
	# For each row i, we earlier defined M, which compares one rare read to all common reads.
	# Now we want to have the following for each compared read at row i:
    #   B[i] = rightmost prefix position such that M[i, 1:B[i]] are all true (0 if M[i,1] is false)
    #   E[i] = leftmost suffix position such that M[i, E[i]:end] are all true (cols+1 if M[i,end] is false)

	rows, cols = size(M)

	B = zeros(Int, rows)

	E = fill(cols + 1, rows)

	# iterate over each row of compared read i
	for i ∈ 1:rows
		# iterate from start to end of unique read i, 
		# and save position j for which all positions left to it (including it)
		# are identical to the rare unique read (i.e. M[i, 1:j] are all true)
		for j ∈ 1:cols
			# keep updating B[i] to be the rightmost possible position,
			# or break the loop if we find a position where M[i, j] is false
			if M[i, j]
				B[i] = j
			else
				break
			end
		end
		# iterate from end to start of compared unique read i,
		# and save position j for which all positions right to it (including it)
		# are identical to the rare unique read (i.e. M[i, j:end] are all true)
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


function are_there_chimeric_reads(B, E)
	return maximum(B) >= minimum(E)
end

function are_there_chimeric_reads(M)
	B, E = find_B_and_E_for_M(M)
	return are_there_chimeric_reads(B, E)
end



function find_all_chimeric_reads(B, E)
	# For each row i, we earlier defined using M (which compares one rare read to all common reads):
    #   B[i] = rightmost prefix position such that M[i, 1:B[i]] are all true (0 if M[i,1] is false)
    #   E[i] = leftmost suffix position such that M[i, E[i]:end] are all true (cols+1 if M[i,end] is false)
    #
    # A "chimeric pair" (r_common, r_common2) exists if E[r_common2] ≤ B[r_common].
    # This implementation returns all such index pairs (original row indices of M).

	# Sort B and E, keep permutations so we can report original row indices
    idxs_B = sortperm(B)
    idxs_E = sortperm(E)
    sorted_B = B[idxs_B]
    sorted_E = E[idxs_E]

	# Two-pointer sweep:
    # As i increases (sorted_B nondecreasing), the set {E <= B[i]} only grows,
    # so we can advance j monotonically.
	pairs = Vector{Tuple{Int, Int}}()
	j = 1
	N = length(sorted_B) # same as length(sorted_E)

	@inbounds for i ∈ 1:N
		Bᵢ = sorted_B[i]
		while j <= N && sorted_E[j] <= Bᵢ
			j += 1
		end
		
		# All E indices in 1:(j-1) satisfy E <= B[i]
		Bᵢ_orig = idxs_B[i]
		for k ∈ 1:j-1
			# println("k: $k")
			push!(pairs, (Bᵢ_orig, idxs_E[k]))
		end
	end

	return pairs
end


function find_all_chimeric_reads(M)
	# For each row i, define:
    #   B[i] = rightmost prefix position such that M[i, 1:B[i]] are all true (0 if M[i,1] is false)
    #   E[i] = leftmost suffix position such that M[i, E[i]:end] are all true (cols+1 if M[i,end] is false)
    #
    # A "chimeric pair" (r_common, r_common2) exists if E[r_common2] ≤ B[r_common].
    # This implementation returns all such index pairs (original row indices of M).

	B, E = find_B_and_E_for_M(M)

	return find_all_chimeric_reads(B, E)
end


extract_chimeric_reads_sites_intersection(all_chimeric_reads_indices, B, E) = [(B[i], E[j]) for (i, j) ∈ all_chimeric_reads_indices]



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





function prepare_expression_df(
	expression_file, 
	x_common_proteins = 300, 
	y_rare_proteins = 1000, 
	discard_reassigned_reads::Bool = true, 
	x_and_y_proteins_denote_fractions::Bool = false
)

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

	num_of_rows = nrow(expression_df)
	
	if x_and_y_proteins_denote_fractions

		if y_rare_proteins + x_common_proteins > 1
			throw(BoundsError("when using `x_and_y_proteins_denote_fractions`, `y_rare_proteins + x_common_proteins` must be less than or equal to 1"))
		end

		# convert to ints, use them instead of `x_common_proteins` and `y_rare_proteins` and also return them to `process_one_sample`
		# for further use downstream
		actual_x_common_proteins = convert(Int, round(x_common_proteins * num_of_rows))
		actual_y_rare_proteins = convert(Int, round(y_rare_proteins * num_of_rows))

		if actual_x_common_proteins == 0 || actual_y_rare_proteins == 0
		throw(
			"converting the x_common_proteins and y_rare_proteins fractions into integers failed - they must be greater than 0 "
			* "(expression_file: $expression_file "
			* "x_common_proteins: $x_common_proteins, y_rare_proteins: $y_rare_proteins "
			* "num_of_rows: $num_of_rows "
			* "actual_x_common_proteins: $actual_x_common_proteins, actual_y_rare_proteins: $actual_y_rare_proteins)"
		)
	end
	
	else

		if y_rare_proteins + x_common_proteins > num_of_rows
			throw(BoundsError("y_rare_proteins + x_common_proteins must be less than or equal to the number of rows in the expression DataFrame"))
		end

		# for compatibility with the previous version of the code, if `x_and_y_proteins_denote_fractions` is false, we just return the input values of `x_common_proteins` and `y_rare_proteins` to `process_one_sample` for further use downstream
		actual_x_common_proteins = x_common_proteins
		actual_y_rare_proteins = y_rare_proteins

	end

	# if actual_x_common_proteins == 0 || actual_y_rare_proteins == 0
	# 	throw(BoundsError("x_common_proteins and y_rare_proteins must be greater than 0"))
	# end
	
	# keep only the y_rare_proteins rarest and x_common_proteins most common distinct proteins
	# (the first y_rare_proteins rows and the last x_common_proteins rows)
	# rare_expression_df = expression_df[1:y_rare_proteins, :]
	# common_expression_df = expression_df[end-x_common_proteins+1:end, :]
	rare_expression_df = expression_df[1:actual_y_rare_proteins, :]
	common_expression_df = expression_df[end-actual_x_common_proteins+1:end, :]



	insertcols!(rare_expression_df, "ExpressionStatus" => "Rare")
	insertcols!(common_expression_df, "ExpressionStatus" => "Common")
	expression_df = vcat(rare_expression_df, common_expression_df)



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

	return expression_df, actual_x_common_proteins, actual_y_rare_proteins
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
	# select!(unique_reads_df, vcat(["Gene", "Transcript"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	select!(unique_reads_df, vcat(["Gene", "Transcript", "NumOfReads"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	rename!(unique_reads_df, "Transcript" => "UniqueRead")
	return unique_reads_df
end


struct ReassignmentMetadata
	Protein::Any
	OriginalOrAdditionalSupportingProtein::Any
	TotalWeightedSupportingReads::Any
end



function process_one_sample(
	platform,
	sample_name,
	expression_file,
	unique_proteins_file, unique_proteins_first_col_pos,
	unique_reads_file, unique_reads_first_col_pos,
	x_common_proteins = 300, y_rare_proteins = 1000,
	discard_reassigned_reads::Bool = true,
	soft_comparison::Bool = false,
	x_and_y_proteins_denote_fractions::Bool = false,
)

	expression_df, actual_x_common_proteins, actual_y_rare_proteins = prepare_expression_df(
		expression_file, x_common_proteins, y_rare_proteins, discard_reassigned_reads, x_and_y_proteins_denote_fractions
	)
	unique_proteins_df = prepare_unique_proteins_df(unique_proteins_file, unique_proteins_first_col_pos)
	unique_reads_df = prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)

	unique_reads_statuses_df = unique_reads_df[:, 4:end]
	nrow(unique(unique_reads_statuses_df)) == nrow(unique_reads_statuses_df) || error(
		"There are duplicate rows in the unique reads editing statuses, which should not happen. Please check the input file: $unique_reads_file"
	)
	editing_sites = parse.(Int, names(unique_reads_statuses_df))
	editing_percents = [
		100 * count(==(1), col) / count(!=(-1), col)
		for col ∈ eachcol(unique_reads_statuses_df)
	]

	unique_reads_and_proteins_df = innerjoin(unique_proteins_df, unique_reads_df, on = ["Gene", "UniqueRead"])

	expression_df = leftjoin(expression_df, unique_reads_and_proteins_df, on = ["Gene", "OriginalOrAdditionalSupportingProtein" => "Protein"])

	# now we have a dataframe with the original distinct proteins (300 most common and 1000 rarest) and their expression levels,
	# as well as the indistinguishable proteins supporting them,
	# and the unique reads supporting them, in turn
	# (each row is a single unique read supporting a unique protein, wether it is a distinct or indistinguishable one which underwent reasignment)

	common_expression_df = expression_df[expression_df.ExpressionStatus.=="Common", :]
	rare_expression_df = expression_df[expression_df.ExpressionStatus.=="Rare", :]

	total_common_unique_reads = nrow(common_expression_df)
	total_rare_unique_reads = nrow(rare_expression_df)
	total_common_reads = sum(common_expression_df.NumOfReads)
	total_rare_reads = sum(rare_expression_df.NumOfReads)

	# new_first_reads_col_pos = 7
	new_first_reads_col_pos = 8  # since adding the number of reads comprising each unique read

	# common_expression_df[:, begin:new_first_reads_col_pos-1]

	# gdf = groupby(common_expression_df[:, ["Gene", "ExpressionStatus", "UniqueRead", "Protein", "OriginalOrAdditionalSupportingProtein", "TotalWeightedSupportingReads"]], ["Gene", "ExpressionStatus", "UniqueRead"])
	gdf = groupby(
		common_expression_df[
			:,
			["Gene", "ExpressionStatus", "UniqueRead", "Protein", "OriginalOrAdditionalSupportingProtein", "TotalWeightedSupportingReads"]
		],
		["Gene", "ExpressionStatus", "UniqueRead"]
	)

	# unique(common_expression_df[!, "UniqueRead"])
	# length(gdf)

	compact_common_expression_df = combine(gdf) do subdf
		(; Metadata = hcat(collect(ReassignmentMetadata.(subdf.Protein, subdf.OriginalOrAdditionalSupportingProtein, subdf.TotalWeightedSupportingReads))))
	end

	# compact_common_expression_df[!, "Metadata"][1]
	# compact_common_expression_df[!, "Metadata"][1][1]

	compact_common_expression_df = leftjoin(
		compact_common_expression_df,
		# unique(common_expression_df[!, new_first_reads_col_pos-1:end], "UniqueRead"),
		unique(common_expression_df[!, new_first_reads_col_pos-2:end], "UniqueRead"), # -2 instead of -1 since we added the "NumOfReads" col
		on = ["UniqueRead"],
	)

	# compact_common_expression_first_reads_col_pos = 5
	compact_common_expression_first_reads_col_pos = 6 # since adding the number of reads comprising each unique read

	# common_unique_reads_editing_status = compact_common_expression_df[:, compact_common_expression_first_reads_col_pos:end]

	results = []

	# iterate over each rare protein and its unique reads, 
	# and compare their editing status to those of the common proteins' unique reads, 
	# to find chimeric reads

	original_rare_proteins = unique(rare_expression_df[!, :Protein])

	# TODO uncomment - this is for testing
	# one_original_rare_protein = original_rare_proteins[1]

	for one_original_rare_protein ∈ original_rare_proteins

		unique_reads_of_one_original_rare_protein_expression_df = rare_expression_df[rare_expression_df.Protein.==one_original_rare_protein, :]
		unique_reads_of_one_original_rare_protein_expression = unique(unique_reads_of_one_original_rare_protein_expression_df[!, "UniqueRead"])

		# TODO uncomment - this is for testing
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

			B, E = find_B_and_E_for_M(M)
			is_chimeric = are_there_chimeric_reads(B, E)
			all_chimeric_reads_indices = find_all_chimeric_reads(B, E)
			chimeric_reads_sites_intersection = extract_chimeric_reads_sites_intersection(all_chimeric_reads_indices, B, E)

			result = (
				Platform = platform,
				Sample = sample_name,
				EditingSites = editing_sites,
				EditingPercents = editing_percents,
				IsSoftComparison = soft_comparison,
				XCommonProteins = x_common_proteins,
				YRareProteins = y_rare_proteins,
				XYProteinsDenoteFractions = x_and_y_proteins_denote_fractions,
				ActualXCommonProteins = actual_x_common_proteins,
				ActualYRareProteins = actual_y_rare_proteins,
				TotalCommonUniqueReads = total_common_unique_reads,
				TotalRareUniqueReads = total_rare_unique_reads,
				TotalCommonReads = total_common_reads,
				TotalRareReads = total_rare_reads,
				ReassignendReadsDiscarded = discard_reassigned_reads,
				Protein = one_original_rare_protein,
				UniqueRead = one_unique_read_of_one_original_rare_protein,
				IsChimeric = is_chimeric,
				ChimericReadsIndices = all_chimeric_reads_indices,
				NumOfChimericCombinations = length(all_chimeric_reads_indices),
				ChimericReadsIntersectingSitesIndices = chimeric_reads_sites_intersection,
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
	platforms,
	samples,
	unique_reads_files,
	unique_proteins_files,
	expression_files,
	unique_reads_first_col_pos,
	unique_proteins_first_col_pos,
	discard_reassigned_reads::Bool = true,
	x_and_y_proteins_denote_fractions::Bool = false,
)
	x_and_y_common_and_rare_proteins = [(x, y) for x in X_common_proteins for y in Y_rare_proteins]

	results_dfs = []

	for soft_comparison in [true, false]

		for (x_common_proteins, y_rare_proteins) ∈ x_and_y_common_and_rare_proteins

			for (platform, sample_name, expression_file, unique_reads_file, unique_proteins_file) ∈ zip(
				platforms, samples, expression_files, unique_reads_files, unique_proteins_files)

				if x_and_y_proteins_denote_fractions
					println(
						"Processing $sample_name ($platform) with soft_comparison = $soft_comparison, "
						* "using $(x_common_proteins * 100)% common proteins and $(y_rare_proteins * 100)% rare proteins..."
					)
				else
					println(
						"Processing $sample_name ($platform) with soft_comparison = $soft_comparison, "
						* "using $x_common_proteins common proteins and $y_rare_proteins rare proteins..."
					)
				end
				
				results_df = try
                    process_one_sample(
						platform,
                        sample_name,
                        expression_file,
                        unique_proteins_file, unique_proteins_first_col_pos,
                        unique_reads_file, unique_reads_first_col_pos,
                        x_common_proteins, y_rare_proteins,
                        discard_reassigned_reads,
                        soft_comparison,
						x_and_y_proteins_denote_fractions
                    )
                catch e
                    if e isa BoundsError
                        @warn "Skipping sample due to BoundsError" platform sample_name soft_comparison x_common_proteins y_rare_proteins exception = (e, catch_backtrace())
                        continue
                    end
                    rethrow()
                end

				chimeric_results_df = results_df[results_df.IsChimeric.==true, :]
				println("Number of chimeric reads in $sample_name: $(nrow(chimeric_results_df))")
				# println(chimeric_results_df)
				push!(results_dfs, results_df)
			end
		end

	end

	results_df = vcat(results_dfs...)

	stats_df = combine(
		groupby(
			results_df,
			[:Platform, :IsSoftComparison, :XYProteinsDenoteFractions, :XCommonProteins, :YRareProteins, :Sample],
		),
		:ActualXCommonProteins => first => :ActualXCommonProteins,
		:ActualYRareProteins => first => :ActualYRareProteins,
		:EditingSites => unique => :EditingSites,
		:EditingPercents => unique => :EditingPercents,
		:IsChimeric => sum => :NumOfChimericReads,
		:TotalCommonUniqueReads => first => :TotalCommonUniqueReads,
		:TotalRareUniqueReads => first => :TotalRareUniqueReads,
		:TotalCommonReads => first => :TotalCommonReads,
		:TotalRareReads => first => :TotalRareReads,
	)
	insertcols!(
		stats_df, 
		"NumOfChimericReads", 
		# "%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.YRareProteins, 
		"%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.ActualYRareProteins, 
		after=true
	)

	return stats_df, results_df
end



function threaded_per_platform_stats_df(
    X_common_proteins,
    Y_rare_proteins,
    platforms,
    samples,
    unique_reads_files,
    unique_proteins_files,
    expression_files,
    unique_reads_first_col_pos,
    unique_proteins_first_col_pos,
    discard_reassigned_reads::Bool = true,
    x_and_y_proteins_denote_fractions::Bool = false,
)
    x_and_y_common_and_rare_proteins = [(x, y) for x in X_common_proteins for y in Y_rare_proteins]

    # Build a flat job list so we can thread over it safely
    jobs = Tuple[]
    for soft_comparison in (true, false)
        for (x_common_proteins, y_rare_proteins) in x_and_y_common_and_rare_proteins
            for (platform, sample_name, expression_file, unique_reads_file, unique_proteins_file) in zip(
                platforms, samples, expression_files, unique_reads_files, unique_proteins_files
            )
                push!(jobs, (
                    soft_comparison,
                    x_common_proteins,
                    y_rare_proteins,
                    platform,
                    sample_name,
                    expression_file,
                    unique_reads_file,
                    unique_proteins_file,
                ))
            end
        end
    end

    out = Vector{Union{Nothing, DataFrame}}(undef, length(jobs))

    Threads.@threads for idx in eachindex(jobs)
        soft_comparison,
        x_common_proteins,
        y_rare_proteins,
        platform,
        sample_name,
        expression_file,
        unique_reads_file,
        unique_proteins_file = jobs[idx]

        # (Optional) avoid println here; it will interleave across threads.

        results_df = try
            process_one_sample(
                platform,
                sample_name,
                expression_file,
                unique_proteins_file, unique_proteins_first_col_pos,
                unique_reads_file, unique_reads_first_col_pos,
                x_common_proteins, y_rare_proteins,
                discard_reassigned_reads,
                soft_comparison,
                x_and_y_proteins_denote_fractions,
            )
        catch e
            if e isa BoundsError
                # skip this job
                out[idx] = nothing
                continue
            end
            rethrow()
        end

        # process_one_sample already includes :Platform and :IsSoftComparison in results_df
        out[idx] = results_df
    end

    results_dfs = (df for df in out if df !== nothing)
    results_df = isempty(results_dfs) ? DataFrame() : vcat(results_dfs...)

    # If everything was skipped, return empty frames with a minimal schema
    if nrow(results_df) == 0
        return DataFrame(), results_df
    end

    stats_df = combine(
        groupby(
            results_df,
            [:Platform, :IsSoftComparison, :XYProteinsDenoteFractions, :XCommonProteins, :YRareProteins, :Sample],
        ),
        :ActualXCommonProteins => first => :ActualXCommonProteins,
        :ActualYRareProteins => first => :ActualYRareProteins,
        :EditingSites => unique => :EditingSites,
        :EditingPercents => unique => :EditingPercents,
        :IsChimeric => sum => :NumOfChimericReads,
        :TotalCommonUniqueReads => first => :TotalCommonUniqueReads,
        :TotalRareUniqueReads => first => :TotalRareUniqueReads,
        :TotalCommonReads => first => :TotalCommonReads,
        :TotalRareReads => first => :TotalRareReads,
    )
    insertcols!(
        stats_df,
        "NumOfChimericReads",
        "%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.ActualYRareProteins,
        after = true,
    )

    return stats_df, results_df
end






pacbio_platforms = vcat(
	fill("PacBio1", 2),
	fill("PacBio2", 2),
	fill("PacBio3", 3),
)
pacbio_samples = [
	"GRIA2", "PCLO",
	"ADAR1", "IQEC1",
	"GRIA2", "ADAR1", "IQEC1",
]
pacbio_unique_reads_files = [
	# PacBio1
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	# PacBio2
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_reads.csv.gz",
	# PacBio3
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.unique_reads.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.unique_reads.csv.gz",
]
pacbio_unique_proteins_files = [
	# PacBio1
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	# PacBio2
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC.Merged.r64296e203404D01.aligned.sorted.MinRQ998.unique_proteins.csv.gz",
	# PacBio3
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141693_c0_seq1.merged.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.unique_proteins.csv.gz",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp141565_c6_seq3.merged.MinRQ998.unique_proteins.csv.gz",
]
pacbio_expression_files = [
	# PacBio1
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio2
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/UMILongReads.MergedSamples/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	# PacBio3
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/GRIA2.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/ADAR1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
	"D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/IQEC1.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv",
]


illumina_samples = [
	"RUSC2",
	"TRIM2",
	"CA2D3",
	"ABL",
	"DGLA",
	"K0513",
	"KCNAS",
	"ACHA4",
	"ANR17",
	"TWK7",
	"SCN1",
	"CACB2",
	"RIMS2",
	"PCLO",
	"DOP1",
	"IQEC1",
	"CSKI1",
	"MTUS2",
	"ROBO2"
]
illumina_chroms = [
    "comp141881_c0_seq3",
    "comp141044_c0_seq2",
    "comp140439_c0_seq1",
    "comp126362_c0_seq1",
    "comp141517_c0_seq1",
    "comp141840_c0_seq2",
    "comp141640_c0_seq1",
    "comp140987_c3_seq1",
    "comp140910_c2_seq1",
    "comp136058_c0_seq1",
    "comp141378_c0_seq7",
    "comp141158_c1_seq2",
    "comp140712_c0_seq3",
    "comp141882_c0_seq14",
    "comp141880_c1_seq3",
    "comp141565_c6_seq3",
    "comp141684_c0_seq1",
    "comp141532_c3_seq11",
    "comp141574_c0_seq3",
]
illumina_unique_reads_files = [
	"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.$chrom.unique_reads.csv"
    for chrom in illumina_chroms
]
illumina_unique_proteins_files = [
	"/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.$chrom.unique_proteins.csv"
    for chrom in illumina_chroms
]
illumina_expression_files = [
	(
		"/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/FinalFixedExpressionFromCloud/" 
		* sample * "/" * sample * ".DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv"
	)
	for sample in illumina_samples
]
illumina_platforms = fill("Illumina1", length(illumina_samples))


platforms = vcat(pacbio_platforms, illumina_platforms)
samples = vcat(pacbio_samples, illumina_samples)
unique_reads_files = vcat(pacbio_unique_reads_files, illumina_unique_reads_files)
unique_proteins_files = vcat(pacbio_unique_proteins_files, illumina_unique_proteins_files)
expression_files = vcat(pacbio_expression_files, illumina_expression_files)


unique_reads_first_col_pos = 9
unique_proteins_first_col_pos = 15

discard_reassigned_reads = true

# X_common_proteins = [300, 500, 1_000]
# Y_rare_proteins = [1_000, 5_000, 10_000]
# x_and_y_proteins_denote_fractions = false

X_common_proteins = [0.01, 0.05, 0.1]
Y_rare_proteins = [0.1, 0.25, 0.5]
x_and_y_proteins_denote_fractions = true



stats_df, results_df = threaded_per_platform_stats_df(
	X_common_proteins,
	Y_rare_proteins,
	platforms,
	samples,
	unique_reads_files,
	unique_proteins_files,
	expression_files,
	unique_reads_first_col_pos,
	unique_proteins_first_col_pos,
	discard_reassigned_reads,
	x_and_y_proteins_denote_fractions
)










# length(samples) * length(X_common_proteins) * length(Y_rare_proteins) * 2 # 2 for soft and strict comparison
# nrow(pacbio_stats_df)


# # Expected number of groups (one row per Platform/Sample/X/Y/Soft)
# expected = length(samples) * length(X_common_proteins) * length(Y_rare_proteins) * 2

# # Actual number of unique groups present
# actual = nrow(unique(pacbio_stats_df, [:Platform, :Sample, :XCommonProteins, :YRareProteins, :IsSoftComparison]))

# println("expected groups = ", expected)
# println("actual groups   = ", actual)
# println("raw nrow(df)    = ", nrow(pacbio_stats_df))


# # Build expected key table
# expected_keys = DataFrame(
#     Platform = platforms,
#     Sample = samples,
# )
# expected_keys = DataFrames.crossjoin(
#     expected_keys,
#     DataFrame(XCommonProteins = X_common_proteins),
#     DataFrame(YRareProteins = Y_rare_proteins),
#     DataFrame(IsSoftComparison = [true, false]),
# )

# have_keys = unique(pacbio_stats_df, [:Platform, :Sample, :XCommonProteins, :YRareProteins, :IsSoftComparison])

# missing_keys = antijoin(expected_keys, have_keys, on = [:Platform, :Sample, :XCommonProteins, :YRareProteins, :IsSoftComparison])

# println("missing rows = ", nrow(missing_keys))
# show(missing_keys, allrows = true, allcols = true)





# # # TODO uncomment to test one sample
# i = 1
# # i = 5
# sample_name = samples[i]
# expression_file = expression_files[i]
# unique_reads_file = unique_reads_files[i]
# unique_proteins_file = unique_proteins_files[i]
# # x_common_proteins = 300
# # y_rare_proteins = 1000
# discard_reassigned_reads = true
# # soft_comparison = false
# soft_comparison = true
# # x_common_proteins = 0.5
# # y_rare_proteins = 0.1
# x_common_proteins = 0.01
# y_rare_proteins = 0.1
# x_and_y_proteins_denote_fractions = true

# results_df = process_one_sample(
#     sample_name,    
#     expression_file, 
#     unique_proteins_file, unique_proteins_first_col_pos,
#     unique_reads_file, unique_reads_first_col_pos,
# 	x_common_proteins, y_rare_proteins,
# 	discard_reassigned_reads,
# 	soft_comparison,
# 	x_and_y_proteins_denote_fractions,
# )





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
	# log2
)



function plot_stats_df_3(
    X_common_proteins,
    Y_rare_proteins,
    stats_df,
    platforms,
    samples,
	x_and_y_proteins_denote_fractions,
    yscale = identity,
    colors = Makie.wong_colors(),
)
    # One bar per (platform,sample) occurrence, in the given order.
    # Color by platform; legend = platforms.
    # Per subplot: only show ticks for samples that exist in that subplot's subdf.
    # Y-axis is synchronized across all subplots within each (soft/strict) figure.
    #
    # Also draws vertical separator lines at *global* platform boundaries that fall within the
    # subplot's shown x-range (even if one side of the boundary has no bars in that subplot).

    # Keep platform order as it appears in `platforms`
    uniq_platforms = String[]
    for p in platforms
        if !(p in uniq_platforms)
            push!(uniq_platforms, p)
        end
    end
    color_by_platform = Dict(p => colors[i] for (i, p) in enumerate(uniq_platforms))

    nrows = length(X_common_proteins)
    ncols = length(Y_rare_proteins)

    # Pre-build global positions for expected (Platform, Sample) occurrences
    positions = Dict{Tuple{Any, Any}, Vector{Int}}()
    for (idx, (p, s)) in enumerate(zip(platforms, samples))
        key = (p, s)
        if haskey(positions, key)
            push!(positions[key], idx)
        else
            positions[key] = [idx]
        end
    end

    for soft_comparison in (true, false)
        fig = Figure(size = (1050, 650))
        axes = Axis[]

        # # y tick step: 20 for soft, 1 for strict
        # ydtick = soft_comparison ? 20 : 1

        # # Precompute a single yupper for the whole figure (this soft/strict)
        # df_sc = stats_df[stats_df.IsSoftComparison .== soft_comparison, :]
        # global_ymax = nrow(df_sc) == 0 ? 0.0 : maximum(df_sc[!, "%OfChimericReads"])
        # global_yupper = global_ymax == 0 ? (soft_comparison ? 100.0 : 5.0) : ceil(global_ymax / ydtick) * ydtick
        # global_yticks = 0:ydtick:global_yupper

        for (i, x_common_proteins) in enumerate(X_common_proteins)
            for (j, y_rare_proteins) in enumerate(Y_rare_proteins)
				if x_and_y_proteins_denote_fractions
					subtitle = "Common: $(format(x_common_proteins * 100, commas=true))%\nRare: $(format(y_rare_proteins * 100, commas=true))%"
				else
                	subtitle = "Common: $(format(x_common_proteins, commas=true))\nRare: $(format(y_rare_proteins, commas=true))"
				end

                subdf = stats_df[
                    stats_df.XCommonProteins .== x_common_proteins .&&
                    stats_df.YRareProteins .== y_rare_proteins .&&
                    stats_df.IsSoftComparison .== soft_comparison,
                    :,
                ]

                # Map each row in this subplot to its global x position (based on platforms/samples order)
                used = Dict{Tuple{Any, Any}, Int}()
                x_for_row = similar(subdf.Platform, Int)
                for r in 1:nrow(subdf)
                    key = (subdf.Platform[r], subdf.Sample[r])
                    k = get!(used, key, 0) + 1
                    used[key] = k
                    if !haskey(positions, key) || k > length(positions[key])
                        x_for_row[r] = length(samples) # fallback
                    else
                        x_for_row[r] = positions[key][k]
                    end
                end

                # Sort bars by x
                perm = sortperm(x_for_row)
                xs = x_for_row[perm]
                ys = subdf[perm, "%OfChimericReads"]
                bar_colors = [color_by_platform[p] for p in subdf[perm, :Platform]]

                # Per-subplot xticks: only positions that actually exist in this subplot
                xticks_sub = if isempty(xs)
                    (Int[], String[])
                else
                    uniq_xs = unique(xs)
                    sort!(uniq_xs)
                    (uniq_xs, [samples[k] for k in uniq_xs])
                end

                ax = Axis(
                    fig[i, j],
                    xticks = xticks_sub,
                    xticklabelrotation = π / 2,
                    subtitle = subtitle,
                    yscale = yscale,
                    xgridvisible = false,
                    xminorgridvisible = false,
                    # yticks = global_yticks,
                )
                push!(axes, ax)

                # # Sync y-limits across all subplots in this figure
                # ylims!(ax, 0, global_yupper)

                barplot!(
                    ax,
                    xs,
                    ys;
                    color = bar_colors,
                    gap = 0.25,
                    strokecolor = :black,
                    strokewidth = 1,
                )

                # Draw separators at global platform boundaries that fall within the shown x-range
                if !isempty(xticks_sub[1])
                    xmin = minimum(xticks_sub[1])
                    xmax = maximum(xticks_sub[1])

                    for k in 1:length(samples)-1
                        if platforms[k] != platforms[k + 1]
                            xsep = k + 0.5
                            if xmin <= xsep <= xmax
                                vlines!(ax, [xsep]; color = (:black, 0.15), linewidth = 1)
                            end
                        end
                    end
                end
            end
        end

        # Legend column (right side): platforms
        legend_col = ncols + 1
        leg_gl = GridLayout()
        fig[1:nrows, legend_col] = leg_gl

        # Tiny hidden axis to host scatter handles
        legax = Axis(leg_gl[1, 1])
        hidedecorations!(legax)
        hidespines!(legax)
        rowsize!(leg_gl, 1, 1)

        platform_handles = [
            scatter!(legax, [0.0], [0.0]; color = color_by_platform[p], markersize = 12)
            for p in uniq_platforms
        ]

        if nrows >= 2
            Legend(leg_gl[2:nrows, 1], platform_handles, uniq_platforms; title = "Platform", framevisible = false)
        else
            Legend(leg_gl[1, 1], platform_handles, uniq_platforms; title = "Platform", framevisible = false)
        end

        linkxaxes!(axes...)
        linkyaxes!(axes...)
        Label(fig[begin:end, 0], "% of chimeric reads", rotation = pi / 2)

        Label(
            fig[0, begin+1:end],
            soft_comparison ? "Soft comparison" : "Strict comparison",
            fontsize = 18,
        )

        display(fig)
    end
end

# ...existing code...




plot_stats_df_3(
	X_common_proteins,
	Y_rare_proteins,
	pacbio_stats_df,
	platforms,
	samples,
	x_and_y_proteins_denote_fractions
	# log2
)


plot_stats_df_3(
	X_common_proteins,
	Y_rare_proteins,
	illumina_stats_df,
	illumina_platforms,
	illumina_samples,
	x_and_y_proteins_denote_fractions
	# log2
)

