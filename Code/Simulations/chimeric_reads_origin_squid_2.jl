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

extract_chimerizing_pairs(common_unique_proteins, all_chimeric_prots_indices) = [
	(common_unique_proteins[i], common_unique_proteins[j])
	for (i, j) ∈ all_chimeric_prots_indices
]


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




function validate_chimerizing_pairs_are_made_of_different_elements(input_file, chimerizing_pairs)
	for (i, (r1, r2)) in enumerate(chimerizing_pairs)
		r1 == r2 && error(
			"Found a chimeric pair which is identical to one another, which should not happen. 
			Please check the input file: $input_file. The first such pair is: ($r1, $r2) (#$i pair)"
		)
	end
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


# function prepare_expression_df_old(
# 	expression_file, 
# 	X_common_proteins, 
# 	Y_rare_proteins,
# 	x_and_y_proteins_denote_fractions::Bool = true
# )

# 	# TODO use an input largest solutuon rather than randomly selecting one

# 	expression_df = DataFrame(CSV.File(expression_file, delim = "\t"))

# 	# count the number of distinct proteins per solution
# 	expression_df = transform(groupby(expression_df, "#Solution"), nrow => "DistinctProteins")

# 	# filter the DataFrame to keep only the rows with the maximum number of distinct proteins
# 	expression_df = subset(expression_df, :DistinctProteins => x -> x .== maximum(x))

# 	# if there are multiple solutions with the same maximum number of distinct proteins, keep only one
# 	max_rand_solution = sample(MersenneTwister(1892), unique(expression_df[!, "#Solution"]), 1, replace = false)
# 	expression_df = subset(expression_df, "#Solution" => x -> x .== max_rand_solution)

# 	# sort the DataFrame by the total expression level, from the lowest to the highest
# 	expression_df = sort(expression_df, "TotalWeightedSupportingReads")

# 	max_x_common_proteins = maximum(X_common_proteins)
# 	max_y_rare_proteins = maximum(Y_rare_proteins)
	
# 	num_of_rows = nrow(expression_df)
	
# 	if x_and_y_proteins_denote_fractions

# 		if max_y_rare_proteins + max_x_common_proteins > 1
# 			throw(BoundsError("when using `x_and_y_proteins_denote_fractions`, `max_y_rare_proteins + max_x_common_proteins` must be less than or equal to 1"))
# 		end

# 		# convert to ints, use them instead of `x_common_proteins` and `y_rare_proteins` and also return them to `process_one_sample`
# 		# for further use downstream
# 		actual_X_common_proteins = convert.(Int, round.(X_common_proteins * num_of_rows))
# 		actual_Y_rare_proteins = convert.(Int, round.(Y_rare_proteins * num_of_rows))

# 		if minimum(actual_X_common_proteins) == 0 || minimum(actual_Y_rare_proteins) == 0
# 			throw(
# 				"converting the X_common_proteins and Y_rare_proteins fractions into integers failed - they must be greater than 0 each "
# 				* "(expression_file: $expression_file, "
# 				* "X_common_proteins: $X_common_proteins, Y_rare_proteins: $Y_rare_proteins, "
# 				* "num_of_rows: $num_of_rows, "
# 				* "actual_X_common_proteins: $actual_X_common_proteins, actual_Y_rare_proteins: $actual_Y_rare_proteins)"
# 			)
# 		end
	
# 	else

# 		if max_y_rare_proteins + max_x_common_proteins > num_of_rows
# 			throw(BoundsError("Any y of Y_rare_proteins + x of X_common_proteins must be less than or equal to the number of rows in the expression DataFrame"))
# 		end

# 		# for compatibility with the previous version of the code, if `x_and_y_proteins_denote_fractions` is false, we just return the input values of `x_common_proteins` and `y_rare_proteins` to `process_one_sample` for further use downstream
# 		actual_X_common_proteins = X_common_proteins
# 		actual_Y_rare_proteins = Y_rare_proteins

# 	end

# 	max_actual_x_common_proteins = maximum(actual_X_common_proteins)
# 	max_actual_y_rare_proteins = maximum(actual_Y_rare_proteins)
	
# 	# keep only the y_rare_proteins rarest and x_common_proteins most common distinct proteins
# 	# (the first y_rare_proteins rows and the last x_common_proteins rows)
# 	# rare_expression_df = expression_df[1:y_rare_proteins, :]
# 	# common_expression_df = expression_df[end-x_common_proteins+1:end, :]
# 	common_expression_df = expression_df[end-max_actual_x_common_proteins+1:end, :]
# 	insertcols!(
# 		common_expression_df,
# 		1,
# 		"AscendingExpressionIndex" => 1:max_actual_x_common_proteins
# 	)
# 	# for the subsets of different common proteins, this col represents the index of the smallest subset of the most common proteins that contains each of the common proteins in the common_expression_df, e.g. if we have 5 common proteins and a certain common protein is in the subset of the 3 most common proteins, then its "SmallestXOrYSubset" value will be 3, as it's contained in the subset of the 3 most common proteins, but not in the subset of the 2 most common proteins.
# 	insertcols!(
# 		common_expression_df,
# 		# 2,
# 		"SmallestXOrYSubset" => [actual_X_common_proteins[searchsortedfirst(actual_X_common_proteins, a)] for a in common_expression_df[:, "AscendingExpressionIndex"]]
# 	)
# 	select!(
# 		common_expression_df,
# 		Not("AscendingExpressionIndex")
# 	)
# 	insertcols!(common_expression_df, "ExpressionStatus" => "Common")
	

# 	rare_expression_df = expression_df[1:max_actual_y_rare_proteins, :]
# 	insertcols!(
# 		rare_expression_df,
# 		1,
# 		"DescendingExpressionIndex" => max_actual_y_rare_proteins:-1:1
# 	)
# 	insertcols!(
# 		rare_expression_df,
# 		# 2,
# 		"SmallestXOrYSubset" => [actual_Y_rare_proteins[searchsortedfirst(actual_Y_rare_proteins, a)] for a in rare_expression_df[:, "DescendingExpressionIndex"]]
# 	)
# 	select!(
# 		rare_expression_df,
# 		Not("DescendingExpressionIndex")
# 	)
# 	insertcols!(rare_expression_df, "ExpressionStatus" => "Rare")
	

# 	expression_df = vcat(rare_expression_df, common_expression_df)


# 	# filter out the columns we don't need
# 	select!(
# 		expression_df,
# 		["Gene", "Protein", "#Solution", "TotalWeightedSupportingReads", "ExpressionStatus", "SmallestXOrYSubset"]
# 	)

# 	# return expression_df, actual_X_common_proteins, actual_Y_rare_proteins, max_actual_x_common_proteins, max_actual_y_rare_proteins
# 	return expression_df, actual_X_common_proteins, actual_Y_rare_proteins
# end


# function prepare_expression_df(
# 	expression_file, 
# 	expression_df::DataFrame,
# 	X_common_proteins, 
# 	Y_rare_proteins,
# 	x_and_y_proteins_denote_fractions::Bool = true
# )
# 	# sort the DataFrame by the total expression level, from the lowest to the highest
# 	expression_df = sort(expression_df, "TotalWeightedSupportingReads")

# 	max_x_common_proteins = maximum(X_common_proteins)
# 	max_y_rare_proteins = maximum(Y_rare_proteins)
	
# 	num_of_rows = nrow(expression_df)
	
# 	if x_and_y_proteins_denote_fractions

# 		if max_y_rare_proteins + max_x_common_proteins > 1
# 			throw(BoundsError("when using `x_and_y_proteins_denote_fractions`, `max_y_rare_proteins + max_x_common_proteins` must be less than or equal to 1"))
# 		end

# 		# convert to ints, use them instead of `x_common_proteins` and `y_rare_proteins` and also return them to `process_one_sample`
# 		# for further use downstream
# 		actual_X_common_proteins = convert.(Int, round.(X_common_proteins * num_of_rows))
# 		actual_Y_rare_proteins = convert.(Int, round.(Y_rare_proteins * num_of_rows))

# 		if minimum(actual_X_common_proteins) == 0 || minimum(actual_Y_rare_proteins) == 0
# 			throw(
# 				"converting the X_common_proteins and Y_rare_proteins fractions into integers failed - they must be greater than 0 each "
# 				* "(expression_file: $expression_file, "
# 				* "X_common_proteins: $X_common_proteins, Y_rare_proteins: $Y_rare_proteins, "
# 				* "num_of_rows: $num_of_rows, "
# 				* "actual_X_common_proteins: $actual_X_common_proteins, actual_Y_rare_proteins: $actual_Y_rare_proteins)"
# 			)
# 		end
	
# 	else

# 		if max_y_rare_proteins + max_x_common_proteins > num_of_rows
# 			throw(BoundsError("Any y of Y_rare_proteins + x of X_common_proteins must be less than or equal to the number of rows in the expression DataFrame"))
# 		end

# 		# for compatibility with the previous version of the code, if `x_and_y_proteins_denote_fractions` is false, we just return the input values of `x_common_proteins` and `y_rare_proteins` to `process_one_sample` for further use downstream
# 		actual_X_common_proteins = X_common_proteins
# 		actual_Y_rare_proteins = Y_rare_proteins

# 	end

# 	max_actual_x_common_proteins = maximum(actual_X_common_proteins)
# 	max_actual_y_rare_proteins = maximum(actual_Y_rare_proteins)
	
# 	# keep only the y_rare_proteins rarest and x_common_proteins most common distinct proteins
# 	# (the first y_rare_proteins rows and the last x_common_proteins rows)
# 	# rare_expression_df = expression_df[1:y_rare_proteins, :]
# 	# common_expression_df = expression_df[end-x_common_proteins+1:end, :]
# 	common_expression_df = expression_df[end-max_actual_x_common_proteins+1:end, :]
# 	insertcols!(
# 		common_expression_df,
# 		1,
# 		"AscendingExpressionIndex" => 1:max_actual_x_common_proteins
# 	)
# 	# for the subsets of different common proteins, this col represents the index of the smallest subset of the most common proteins that contains each of the common proteins in the common_expression_df, e.g. if we have 5 common proteins and a certain common protein is in the subset of the 3 most common proteins, then its "SmallestXOrYSubset" value will be 3, as it's contained in the subset of the 3 most common proteins, but not in the subset of the 2 most common proteins.
# 	insertcols!(
# 		common_expression_df,
# 		# 2,
# 		"SmallestXOrYSubset" => [actual_X_common_proteins[searchsortedfirst(actual_X_common_proteins, a)] for a in common_expression_df[:, "AscendingExpressionIndex"]]
# 	)
# 	select!(
# 		common_expression_df,
# 		Not("AscendingExpressionIndex")
# 	)
# 	insertcols!(common_expression_df, "ExpressionStatus" => "Common")
	

# 	rare_expression_df = expression_df[1:max_actual_y_rare_proteins, :]
# 	insertcols!(
# 		rare_expression_df,
# 		1,
# 		"DescendingExpressionIndex" => max_actual_y_rare_proteins:-1:1
# 	)
# 	insertcols!(
# 		rare_expression_df,
# 		# 2,
# 		"SmallestXOrYSubset" => [actual_Y_rare_proteins[searchsortedfirst(actual_Y_rare_proteins, a)] for a in rare_expression_df[:, "DescendingExpressionIndex"]]
# 	)
# 	select!(
# 		rare_expression_df,
# 		Not("DescendingExpressionIndex")
# 	)
# 	insertcols!(rare_expression_df, "ExpressionStatus" => "Rare")
	

# 	expression_df = vcat(rare_expression_df, common_expression_df)


# 	# filter out the columns we don't need
# 	select!(
# 		expression_df,
# 		["Gene", "Protein", "#Solution", "TotalWeightedSupportingReads", "ExpressionStatus", "SmallestXOrYSubset"]
# 	)

# 	# return expression_df, actual_X_common_proteins, actual_Y_rare_proteins, max_actual_x_common_proteins, max_actual_y_rare_proteins
# 	return expression_df, actual_X_common_proteins, actual_Y_rare_proteins
# end


function prepare_expression_df(
	expression_file, 
	expression_df::DataFrame,
	X_common_proteins, 
	Y_rare_proteins,
	x_and_y_proteins_denote_fractions::Bool = true
)
	max_x_common_proteins = maximum(X_common_proteins)
	max_y_rare_proteins = maximum(Y_rare_proteins)
	
	num_of_rows = nrow(expression_df)
	
	if x_and_y_proteins_denote_fractions

		if max_y_rare_proteins + max_x_common_proteins > 1
			throw(BoundsError("when using `x_and_y_proteins_denote_fractions`, `max_y_rare_proteins + max_x_common_proteins` must be less than or equal to 1"))
		end

		# convert to ints, use them instead of `x_common_proteins` and `y_rare_proteins` and also return them to `process_one_sample`
		# for further use downstream
		actual_X_common_proteins = convert.(Int, round.(X_common_proteins * num_of_rows))
		actual_Y_rare_proteins = convert.(Int, round.(Y_rare_proteins * num_of_rows))

		@assert length(actual_X_common_proteins) == length(X_common_proteins)
		@assert length(actual_Y_rare_proteins) == length(Y_rare_proteins)

		if minimum(actual_X_common_proteins) == 0 || minimum(actual_Y_rare_proteins) == 0
			throw(
				"converting the X_common_proteins and Y_rare_proteins fractions into integers failed - they must be greater than 0 each "
				* "(expression_file: $expression_file, "
				* "X_common_proteins: $X_common_proteins, Y_rare_proteins: $Y_rare_proteins, "
				* "num_of_rows: $num_of_rows, "
				* "actual_X_common_proteins: $actual_X_common_proteins, actual_Y_rare_proteins: $actual_Y_rare_proteins)"
			)
		end
	
	else

		if max_y_rare_proteins + max_x_common_proteins > num_of_rows
			throw(BoundsError("Any y of Y_rare_proteins + x of X_common_proteins must be less than or equal to the number of rows in the expression DataFrame"))
		end

		# for compatibility with the previous version of the code, if `x_and_y_proteins_denote_fractions` is false, we just return the input values of `x_common_proteins` and `y_rare_proteins` to `process_one_sample` for further use downstream
		actual_X_common_proteins = X_common_proteins
		actual_Y_rare_proteins = Y_rare_proteins

	end

	max_actual_x_common_proteins = maximum(actual_X_common_proteins)
	max_actual_y_rare_proteins = maximum(actual_Y_rare_proteins)

	# retain only required columns
	select!(
		expression_df,
		["Gene", "Protein", "#Solution", "TotalWeightedSupportingReads",]
	)
	
	# sort the DataFrame by the total expression level, from the lowest to the highest
	expression_df = sort(expression_df, "TotalWeightedSupportingReads")

	insertcols!(
		expression_df,
		"DescendingRareIndex" => 1:num_of_rows,
		"DescendingCommonIndex" => num_of_rows:-1:1
	)

	# expression_df = expression_df[
	# 	(expression_df[!, "DescendingRareIndex"] .<= max_actual_y_rare_proteins
	# 	.|| expression_df[!, "DescendingCommonIndex"] .<= max_actual_x_common_proteins),
	# 	begin:end
	# ]

	expression_df = expression_df[
		(expression_df[!, "DescendingRareIndex"] .<= max_actual_y_rare_proteins) .|
    	(expression_df[!, "DescendingCommonIndex"] .<= max_actual_x_common_proteins),
		begin:end
	]



	# # keep only the y_rare_proteins rarest and x_common_proteins most common distinct proteins
	# # (the first y_rare_proteins rows and the last x_common_proteins rows)
	# # rare_expression_df = expression_df[1:y_rare_proteins, :]
	# # common_expression_df = expression_df[end-x_common_proteins+1:end, :]
	# common_expression_df = expression_df[end-max_actual_x_common_proteins+1:end, :]
	# insertcols!(
	# 	common_expression_df,
	# 	1,
	# 	"AscendingExpressionIndex" => 1:max_actual_x_common_proteins
	# )
	# # for the subsets of different common proteins, this col represents the index of the smallest subset of the most common proteins that contains each of the common proteins in the common_expression_df, e.g. if we have 5 common proteins and a certain common protein is in the subset of the 3 most common proteins, then its "SmallestXOrYSubset" value will be 3, as it's contained in the subset of the 3 most common proteins, but not in the subset of the 2 most common proteins.
	# insertcols!(
	# 	common_expression_df,
	# 	# 2,
	# 	"SmallestXOrYSubset" => [actual_X_common_proteins[searchsortedfirst(actual_X_common_proteins, a)] for a in common_expression_df[:, "AscendingExpressionIndex"]]
	# )
	# select!(
	# 	common_expression_df,
	# 	Not("AscendingExpressionIndex")
	# )
	# insertcols!(common_expression_df, "ExpressionStatus" => "Common")
	

	# rare_expression_df = expression_df[1:max_actual_y_rare_proteins, :]
	# insertcols!(
	# 	rare_expression_df,
	# 	1,
	# 	"DescendingExpressionIndex" => max_actual_y_rare_proteins:-1:1
	# )
	# insertcols!(
	# 	rare_expression_df,
	# 	# 2,
	# 	"SmallestXOrYSubset" => [actual_Y_rare_proteins[searchsortedfirst(actual_Y_rare_proteins, a)] for a in rare_expression_df[:, "DescendingExpressionIndex"]]
	# )
	# select!(
	# 	rare_expression_df,
	# 	Not("DescendingExpressionIndex")
	# )
	# insertcols!(rare_expression_df, "ExpressionStatus" => "Rare")
	

	# expression_df = vcat(rare_expression_df, common_expression_df)


	# # filter out the columns we don't need
	# select!(
	# 	expression_df,
	# 	["Gene", "Protein", "#Solution", "TotalWeightedSupportingReads", "ExpressionStatus", "SmallestXOrYSubset"]
	# )

	# return expression_df, actual_X_common_proteins, actual_Y_rare_proteins, max_actual_x_common_proteins, max_actual_y_rare_proteins
	return expression_df, actual_X_common_proteins, actual_Y_rare_proteins
end


function prepare_expression_df(
	expression_file,
	max_solution::Int, 
	X_common_proteins, 
	Y_rare_proteins,
	x_and_y_proteins_denote_fractions::Bool = true
)
	expression_df = DataFrame(CSV.File(expression_file, delim = "\t"))

	# filter the DataFrame to keep only the rows with the maximum number of distinct proteins, as defined by the input `max_solution` argument
	expression_df = subset(expression_df, "#Solution" => x -> x .== max_solution)

	return prepare_expression_df(
		expression_file, 
		expression_df,
		X_common_proteins, 
		Y_rare_proteins,
		x_and_y_proteins_denote_fractions
	)

end


# function prepare_expression_df(
# 	expression_file, 
# 	X_common_proteins, 
# 	Y_rare_proteins,
# 	x_and_y_proteins_denote_fractions::Bool = true
# )
# 	expression_df = DataFrame(CSV.File(expression_file, delim = "\t"))

# 	# count the number of distinct proteins per solution
# 	expression_df = transform(groupby(expression_df, "#Solution"), nrow => "DistinctProteins")

# 	# filter the DataFrame to keep only the rows with the maximum number of distinct proteins
# 	expression_df = subset(expression_df, :DistinctProteins => x -> x .== maximum(x))

# 	# if there are multiple solutions with the same maximum number of distinct proteins, keep only one
# 	max_rand_solution = sample(MersenneTwister(1892), unique(expression_df[!, "#Solution"]), 1, replace = false)
# 	expression_df = subset(expression_df, "#Solution" => x -> x .== max_rand_solution)

# 	return prepare_expression_df(
# 		expression_file, 
# 		expression_df,
# 		X_common_proteins, 
# 		Y_rare_proteins,
# 		x_and_y_proteins_denote_fractions
# 	)
# end



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
	rename!(unique_reads_df, "Transcript" => "UniqueRead")
	select!(unique_reads_df, vcat(["Gene", "UniqueRead", "NumOfReads"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	new_unique_reads_first_col_pos = 4

	# convert the editing status columns to Int8 to save memory, as they can only take values of 1, 0, or -1
	df_1 = unique_reads_df[:, 1:new_unique_reads_first_col_pos-1]
	df_2 = unique_reads_df[:, new_unique_reads_first_col_pos:end]
	df_2 = convert.(
		Int8,
		df_2
	)
	unique_reads_df = hcat(df_1, df_2)

	return unique_reads_df, new_unique_reads_first_col_pos
end



function compare_editing_statuses(common_status_array, rare_status_array, soft_comparison)
	# Compare the rare read to all common reads (row-wise)
	if soft_comparison
		M = softcomparison.(
			eachrow(common_status_array), 
			Ref(rare_status_array)
		)
		M = hcat(M...)'
	else
		M = common_status_array .== rare_status_array'
	end

	return M
end


function define_out_file(out_dir, platform, sample_name, X_common_proteins::Vector, Y_rare_proteins::Vector)
	return (
		out_dir 
		* "/" 
		* "$(platform).$(sample_name).X_$(join(string.(X_common_proteins), "_")).Y_$(join(string.(Y_rare_proteins), "_")).csv.gz"
	)
end


function save_results_df(
	out_dir,
	results_df,
	platform,
	sample_name,
	X_common_proteins,
	Y_rare_proteins,
)
	results_df = deepcopy(results_df)

	for col in [
		:ChimerizingProteinPairs, 
		:ChimerizingReadPairs, 
		:ChimerizingReadPairsIntersectingSitesIndices
	]
		results_df[!, col] = map(
			v -> join(string.(first.(v), ",", last.(v)), ";"),
			results_df[!, col]
		)
	end

	mkpath(out_dir)  # make absolutely sure it exists (good under threads)
	
	# soft_comparison_interfix_str = soft_comparison ? "soft_comparison" : "strict_comparison"
	# out_file = out_dir * "/" * basename(replace(original_in_file, ".unique_reads.csv.gz" => ".chimeric_reads.$soft_comparison_interfix_str.csv.gz"))
	out_file = define_out_file(out_dir, platform, sample_name, X_common_proteins, Y_rare_proteins)
	
	# Atomic write: write to temp then move into place.
	# tmp_file = out_file * ".tmp." * string(getpid()) * "." * string(threadid())
	tmp_file = try
		out_file * ".tmp." * string(getpid()) * "." * string(threadid())
	catch e
		out_file * ".tmp." * string(getpid()) # if threadid not defined
	end

	# bump buffer to handle very long rows (large joined string)
	# Write to temp first, then move into place (atomic within same filesystem).
    CSV.write(tmp_file, results_df, delim = "\t"; compress = true, bufsize = 64 * 1024 * 1024)
	mv(tmp_file, out_file; force = true)
    
	return out_file
end


"""
Return chimerizing tuples filtered to those where *both* proteins are in `allowed_common_proteins`.
"""
function filter_chimerizing_by_common(
    chimerizing_protein_pairs,
    chimerizing_read_pairs,
    chimerizing_reads_sites_intersection,
    allowed_common_proteins,
)
    # Use Set for O(1) membership
    common_set = allowed_common_proteins isa Set ? allowed_common_proteins : Set(allowed_common_proteins)

    filtered_protein_pairs = eltype(chimerizing_protein_pairs)[]
    filtered_read_pairs = eltype(chimerizing_read_pairs)[]
    filtered_sites_intersection = eltype(chimerizing_reads_sites_intersection)[]

    for (protein_pair, read_pair, sites_intersection) ∈ zip(
        chimerizing_protein_pairs, chimerizing_read_pairs, chimerizing_reads_sites_intersection
    )
        if (protein_pair[1] ∈ common_set) && (protein_pair[2] ∈ common_set)
            push!(filtered_protein_pairs, protein_pair)
            push!(filtered_read_pairs, read_pair)
            push!(filtered_sites_intersection, sites_intersection)
        end
    end

    return filtered_protein_pairs, filtered_read_pairs, filtered_sites_intersection
end


"""
Build a NamedTuple result row (the schema you write to disk).
"""
function make_result_row(;
    platform,
    sample_name,
    max_solution,
    soft_comparison,
    x_common,
    y_rare,
    x_and_y_proteins_denote_fractions,
    actual_x_common,
    actual_y_rare,
    rare_protein,
    one_unique_read_of_rare_protein,
    filtered_protein_pairs,
    filtered_read_pairs,
    filtered_sites_intersection,
)
    num_chim = length(filtered_protein_pairs)
    is_chim = num_chim > 0

    return (
        Platform = platform,
        Sample = sample_name,
        Solution = max_solution,
        IsSoftComparison = soft_comparison,
        XCommonProteins = x_common,
        YRareProteins = y_rare,
        XYProteinsDenoteFractions = x_and_y_proteins_denote_fractions,
        ActualXCommonProteins = actual_x_common,
        ActualYRareProteins = actual_y_rare,
        Protein = rare_protein,
        UniqueRead = one_unique_read_of_rare_protein,
        IsChimeric = is_chim,
        NumOfChimericCombinations = num_chim,
        ChimerizingProteinPairs = filtered_protein_pairs,
        ChimerizingReadPairs = filtered_read_pairs,
        ChimerizingReadPairsIntersectingSitesIndices = filtered_sites_intersection,
    )
end


function process_one_sample(
	platform,
	sample_name,
	expression_file,
	max_solution,
	unique_proteins_file, 
	unique_proteins_first_col_pos,
	unique_reads_file, 
	unique_reads_first_col_pos,
	X_common_proteins, 
	Y_rare_proteins,
	x_and_y_proteins_denote_fractions::Bool = true,
)

	expression_df, actual_X_common_proteins, actual_Y_rare_proteins = prepare_expression_df(
		expression_file, 
		max_solution,
		X_common_proteins, 
		Y_rare_proteins,
		x_and_y_proteins_denote_fractions
	)
	unique_proteins_df = prepare_unique_proteins_df(unique_proteins_file, unique_proteins_first_col_pos)
	unique_reads_df, new_unique_reads_first_col_pos = prepare_unique_reads_df(unique_reads_file, unique_reads_first_col_pos)

	# unique_reads_statuses_df = unique_reads_df[:, new_unique_reads_first_col_pos:end]
	# nrow(unique(unique_reads_statuses_df)) == nrow(unique_reads_statuses_df) || error(
	# 	"There are duplicate rows in the unique reads editing statuses, which should not happen. Please check the input file: $unique_reads_file"
	# )
	# editing_sites = parse.(Int, names(unique_reads_statuses_df))
	# editing_percents = [
	# 	100 * count(==(1), col) / count(!=(-1), col)
	# 	for col ∈ eachcol(unique_reads_statuses_df)
	# ]

	unique_proteins_id_cols = names(unique_proteins_df) # all its cols are id cols
	unique_reads_id_cols = names(unique_reads_df)[1:new_unique_reads_first_col_pos-1]
	num_of_id_cols_in_prots_but_not_in_reads = length(
		setdiff(unique_proteins_id_cols, unique_reads_id_cols)
	)
	new_unique_reads_and_proteins_first_col_pos = (
		new_unique_reads_first_col_pos + num_of_id_cols_in_prots_but_not_in_reads
	)
	
	# a flat df where each row is a protein and one of its unique reads
	unique_reads_and_proteins_df = innerjoin(
		unique_proteins_df, 
		unique_reads_df, 
		on = ["Gene", "UniqueRead"]
	)

	expression_id_cols = names(expression_df)
	unique_reads_and_proteins_id_cols = names(unique_reads_and_proteins_df)[1:new_unique_reads_and_proteins_first_col_pos-1]
	num_of_id_cols_in_prots_and_reads_df_but_not_in_expression_df = length(
		setdiff(unique_reads_and_proteins_id_cols, expression_id_cols)
	)
	new_first_reads_col_pos = (
		size(expression_df, 2) + num_of_id_cols_in_prots_and_reads_df_but_not_in_expression_df + 1
	)
	expression_df = leftjoin(
		expression_df, 
		unique_reads_and_proteins_df, 
		on = ["Gene", "Protein"]
	)

	max_X_common_proteins = maximum(X_common_proteins)
	max_Y_rare_proteins = maximum(Y_rare_proteins)
	max_actual_x_common_proteins = maximum(actual_X_common_proteins)
	max_actual_y_rare_proteins = maximum(actual_Y_rare_proteins)

	X_common_proteins_wo_max_x = setdiff(X_common_proteins, [max_X_common_proteins])
	Y_rare_proteins_wo_max_y = setdiff(Y_rare_proteins, [max_Y_rare_proteins])
	actual_X_common_proteins_wo_max_x = setdiff(actual_X_common_proteins, [max_actual_x_common_proteins])
	actual_Y_rare_proteins_wo_max_y = setdiff(actual_Y_rare_proteins, [max_actual_y_rare_proteins])
	
	# sort in descending order, so that we start with the largest subset of common proteins first, which is more likely to explain the rare proteins and thus save time by skipping the smaller subsets of common proteins that are contained in it
	sort!(X_common_proteins_wo_max_x, rev = true)
	sort!(actual_X_common_proteins_wo_max_x, rev = true)
	# sort in descending order, so that we start with the largest subset of rare proteins first, which is more likely to be explained by the common proteins and thus save time by skipping the smaller subsets of rare proteins that are contained in it
	sort!(Y_rare_proteins_wo_max_y, rev = true)
	sort!(actual_Y_rare_proteins_wo_max_y, rev = true)
	
	# now we have a dataframe with the original distinct proteins (x most common and y rarest) 
	# and their expression levels,
	# and the unique reads supporting them,
	# (each row is a single unique read originally-supporting a unique protein)

	# common_expression_df = expression_df[expression_df.ExpressionStatus.=="Common", :]
	# rare_expression_df = expression_df[expression_df.ExpressionStatus.=="Rare", :]

	common_expression_df = expression_df[expression_df.DescendingCommonIndex .<= max_actual_x_common_proteins, :]
	rare_expression_df = expression_df[expression_df.DescendingRareIndex .<= max_actual_y_rare_proteins, :]

	# total_common_unique_reads = nrow(common_expression_df)
	# total_rare_unique_reads = nrow(rare_expression_df)
	# total_common_reads = sum(common_expression_df.NumOfReads)
	# total_rare_reads = sum(rare_expression_df.NumOfReads)

	common_proteins = common_expression_df[!, :Protein]
	common_unique_reads = common_expression_df[!, :UniqueRead]
	common_unique_reads_editing_status_array = Matrix{Int8}(
		common_expression_df[!, new_first_reads_col_pos:end]
	)

	rare_proteins = rare_expression_df[!, :Protein]
	rare_unique_reads = rare_expression_df[!, :UniqueRead]

	
	# prepare subsets of metadata dfs for the different subsets of common and rare proteins, 
	# to save time by not having to subset them repeatedly inside the loop over the unique reads of the rare proteins
	subset_rare_expression_dfs = Dict(
		y => rare_expression_df[
			rare_expression_df.DescendingRareIndex .<= y, 
			begin:new_first_reads_col_pos-1
		]
		for y in actual_Y_rare_proteins_wo_max_y
	)
	# Before:
    # subset_common_expression_dfs = Dict(
    # 	x => common_expression_df[
    # 		common_expression_df.DescendingCommonIndex .<= x, 
    # 		begin:new_first_reads_col_pos-1
    # 	]
    # 	for x in actual_X_common_proteins_wo_max_x
    # )
    # After: build for *all* X (including max X)
    subset_common_expression_dfs = Dict(
        x => common_expression_df[
            common_expression_df.DescendingCommonIndex .<= x,
            begin:new_first_reads_col_pos-1
        ]
        for x in actual_X_common_proteins
    )
	
	for y in actual_Y_rare_proteins_wo_max_y
		# println("$y, $(nrow(unique(subset_rare_expression_dfs[y], :Protein)))")
		@assert nrow(unique(subset_rare_expression_dfs[y], :Protein)) == y
	end
	# for x in actual_X_common_proteins_wo_max_x
	for x in actual_X_common_proteins
		# println("$x, $(nrow(unique(subset_common_expression_dfs[x], :Protein)))")
		@assert nrow(unique(subset_common_expression_dfs[x], :Protein)) == x
	end

	results = []

	# iterate over each rare protein and its unique reads, 
	# and compare their editing status to those of the common proteins' unique reads, 
	# to find chimeric reads,
	# using both soft and hard comparison, and save the results in a DataFrame

	for (rare_protein, one_unique_read_of_rare_protein) ∈ zip(rare_proteins, rare_unique_reads)

		# rare_protein = rare_proteins[1]
		# one_unique_read_of_rare_protein = rare_unique_reads[1]

		# rare_protein = rare_proteins[2]
		# one_unique_read_of_rare_protein = rare_unique_reads[2]

		# take the one and only row of that unique read
		one_unique_read_of_rare_protein_df = only(
			rare_expression_df[rare_expression_df.UniqueRead.==one_unique_read_of_rare_protein, :]
		)

		one_unique_read_of_rare_protein_editing_status_array = Vector{Int8}(
			one_unique_read_of_rare_protein_df[new_first_reads_col_pos:end]
		)

		for soft_comparison ∈ [true, false]
		
			# soft_comparison = true
			# soft_comparison = false

			# Compare the rare read to all common reads (row-wise)
			M = compare_editing_statuses(
				common_unique_reads_editing_status_array, 
				one_unique_read_of_rare_protein_editing_status_array, 
				soft_comparison
			)

			B, E = find_B_and_E_for_M(M)
			# these are the indices of the max_actual_x_common_proteins that can chimerize to 
			# form the rare protein which is included in the max_actual_y_rare_proteins.
			indices_of_chimeric_read_pairs = find_all_chimeric_reads(B, E)
			chimerizing_protein_pairs = extract_chimerizing_pairs(common_proteins, indices_of_chimeric_read_pairs)
			chimerizing_read_pairs = extract_chimerizing_pairs(common_unique_reads, indices_of_chimeric_read_pairs)
			validate_chimerizing_pairs_are_made_of_different_elements(expression_file, chimerizing_read_pairs)
			chimerizing_reads_sites_intersection = extract_chimeric_reads_sites_intersection(indices_of_chimeric_read_pairs, B, E)
			num_of_chimeric_combinations = length(chimerizing_protein_pairs)
			is_chimeric = num_of_chimeric_combinations > 0

			# this is the result for the max_actual_x_common_proteins and max_actual_y_rare_proteins
			result = (
				Platform = platform,
				Sample = sample_name,
				Solution = max_solution,
				IsSoftComparison = soft_comparison,
				XCommonProteins = max_X_common_proteins,
				YRareProteins = max_Y_rare_proteins,
				XYProteinsDenoteFractions = x_and_y_proteins_denote_fractions,
				ActualXCommonProteins = max_actual_x_common_proteins,
				ActualYRareProteins = max_actual_y_rare_proteins,
				Protein = rare_protein,
				UniqueRead = one_unique_read_of_rare_protein,
				IsChimeric = is_chimeric,
				NumOfChimericCombinations = num_of_chimeric_combinations,
				ChimerizingProteinPairs = chimerizing_protein_pairs,
				ChimerizingReadPairs = chimerizing_read_pairs,
				ChimerizingReadPairsIntersectingSitesIndices = chimerizing_reads_sites_intersection,
			)
			push!(results, result)

			# now we should subset this result for the smaller subsets of common and rare proteins that are contained 
			# in the max_actual_x_common_proteins and max_actual_y_rare_proteins, respectively

			
			# Emit subsets over X when Y is max (currently missing).
            # This makes (x<max_x, y=max_y) exist in the output.
            for (x_common, actual_x_common) ∈ zip(X_common_proteins_wo_max_x, actual_X_common_proteins_wo_max_x)

                subset_common_expression_df = subset_common_expression_dfs[actual_x_common]
                subset_common_proteins = subset_common_expression_df[!, :Protein]

                # subset_chimerizing_protein_pairs = eltype(chimerizing_protein_pairs)[]
                # subset_chimerizing_read_pairs = eltype(chimerizing_read_pairs)[]
                # subset_chimerizing_reads_sites_intersection = eltype(chimerizing_reads_sites_intersection)[]

                # for (protein_pair, read_pair, sites_intersection) ∈ zip(
                #     chimerizing_protein_pairs, chimerizing_read_pairs, chimerizing_reads_sites_intersection
                # )
                #     if (protein_pair[1] ∈ subset_common_proteins) && (protein_pair[2] ∈ subset_common_proteins)
                #         push!(subset_chimerizing_protein_pairs, protein_pair)
                #         push!(subset_chimerizing_read_pairs, read_pair)
                #         push!(subset_chimerizing_reads_sites_intersection, sites_intersection)
                #     end
                # end

                # subset_num_of_chimeric_combinations = length(subset_chimerizing_protein_pairs)
                # subset_is_chimeric = subset_num_of_chimeric_combinations > 0

                # subset_result = (
                #     Platform = platform,
                #     Sample = sample_name,
                #     Solution = max_solution,
                #     IsSoftComparison = soft_comparison,
                #     XCommonProteins = x_common,
                #     YRareProteins = max_Y_rare_proteins,              # <-- key change
                #     XYProteinsDenoteFractions = x_and_y_proteins_denote_fractions,
                #     ActualXCommonProteins = actual_x_common,
                #     ActualYRareProteins = max_actual_y_rare_proteins,  # <-- key change
                #     Protein = rare_protein,
                #     UniqueRead = one_unique_read_of_rare_protein,
                #     IsChimeric = subset_is_chimeric,
                #     NumOfChimericCombinations = subset_num_of_chimeric_combinations,
                #     ChimerizingProteinPairs = subset_chimerizing_protein_pairs,
                #     ChimerizingReadPairs = subset_chimerizing_read_pairs,
                #     ChimerizingReadPairsIntersectingSitesIndices = subset_chimerizing_reads_sites_intersection,
                # )

				fpairs, freads, fsites = filter_chimerizing_by_common(
					chimerizing_protein_pairs,
					chimerizing_read_pairs,
					chimerizing_reads_sites_intersection,
					subset_common_proteins,
				)

				subset_result = make_result_row(
					platform = platform,
					sample_name = sample_name,
					max_solution = max_solution,
					soft_comparison = soft_comparison,
					x_common = x_common,
					y_rare = max_Y_rare_proteins,
					x_and_y_proteins_denote_fractions = x_and_y_proteins_denote_fractions,
					actual_x_common = actual_x_common,
					actual_y_rare = max_actual_y_rare_proteins,
					rare_protein = rare_protein,
					one_unique_read_of_rare_protein = one_unique_read_of_rare_protein,
					filtered_protein_pairs = fpairs,
					filtered_read_pairs = freads,
					filtered_sites_intersection = fsites,
				)

                push!(results, subset_result)
				
            end

			
			# now iterate over the subsets of rare proteins - if this rare protein isn't found in a subset of 
			# rare proteins, we shouldn't check which subset of common proteins can chimerize it
			for (y_rare, actual_y_rare) ∈ zip(Y_rare_proteins_wo_max_y, actual_Y_rare_proteins_wo_max_y)
				
				# y_rare, actual_y_rare = Y_rare_proteins_wo_max_y[1], actual_Y_rare_proteins_wo_max_y[1]
				
				subset_rare_expression_df = subset_rare_expression_dfs[actual_y_rare]
				
				subset_rare_proteins = subset_rare_expression_df[!, :Protein]
				
				if rare_protein ∉ subset_rare_proteins
					continue
				end

				# Before:
				# for (x_common, actual_x_common) ∈ zip(X_common_proteins_wo_max_x, actual_X_common_proteins_wo_max_x)
				
				# After (include x=max too):
				for (x_common, actual_x_common) ∈ zip(X_common_proteins, actual_X_common_proteins)
				
					# x_common, actual_x_common = X_common_proteins_wo_max_x[1], actual_X_common_proteins_wo_max_x[1]
								
					subset_common_expression_df = subset_common_expression_dfs[actual_x_common]

					subset_common_proteins = subset_common_expression_df[!, :Protein]

					# subset_chimerizing_protein_pairs = eltype(chimerizing_protein_pairs)[]
                    # subset_chimerizing_read_pairs = eltype(chimerizing_read_pairs)[]
                    # subset_chimerizing_reads_sites_intersection = eltype(chimerizing_reads_sites_intersection)[]
					
					# for (protein_pair, read_pair, sites_intersection) ∈ zip(
					# 	chimerizing_protein_pairs, chimerizing_read_pairs, chimerizing_reads_sites_intersection
					# )
					# 	if (protein_pair[1] ∈ subset_common_proteins) && (protein_pair[2] ∈ subset_common_proteins)
					# 		push!(subset_chimerizing_protein_pairs, protein_pair)
					# 		push!(subset_chimerizing_read_pairs, read_pair)
					# 		push!(subset_chimerizing_reads_sites_intersection, sites_intersection)
					# 	end
					# end
					
					# subset_num_of_chimeric_combinations = length(subset_chimerizing_protein_pairs)
					# subset_is_chimeric = subset_num_of_chimeric_combinations > 0

					# # this is the result for the current subset of common and rare proteins
					# result = (
					# 	Platform = platform,
					# 	Sample = sample_name,
					# 	Solution = max_solution,
					# 	IsSoftComparison = soft_comparison,
					# 	XCommonProteins = x_common,
					# 	YRareProteins = y_rare,
					# 	XYProteinsDenoteFractions = x_and_y_proteins_denote_fractions,
					# 	ActualXCommonProteins = actual_x_common,
					# 	ActualYRareProteins = actual_y_rare,
					# 	Protein = rare_protein,
					# 	UniqueRead = one_unique_read_of_rare_protein,
					# 	IsChimeric = subset_is_chimeric,
					# 	NumOfChimericCombinations = subset_num_of_chimeric_combinations,
					# 	ChimerizingProteinPairs = subset_chimerizing_protein_pairs,
					# 	ChimerizingReadPairs = subset_chimerizing_read_pairs,
					# 	ChimerizingReadPairsIntersectingSitesIndices = subset_chimerizing_reads_sites_intersection,
					# )

					# push!(results, result)

					fpairs, freads, fsites = filter_chimerizing_by_common(
						chimerizing_protein_pairs,
						chimerizing_read_pairs,
						chimerizing_reads_sites_intersection,
						subset_common_proteins,
					)

					subset_result = make_result_row(
						platform = platform,
						sample_name = sample_name,
						max_solution = max_solution,
						soft_comparison = soft_comparison,
						x_common = x_common,
						y_rare = y_rare,
						x_and_y_proteins_denote_fractions = x_and_y_proteins_denote_fractions,
						actual_x_common = actual_x_common,
						actual_y_rare = actual_y_rare,
						rare_protein = rare_protein,
						one_unique_read_of_rare_protein = one_unique_read_of_rare_protein,
						filtered_protein_pairs = fpairs,
						filtered_read_pairs = freads,
						filtered_sites_intersection = fsites,
					)

					push!(results, subset_result)
				end

			end

		end

	end

	results_df = DataFrame(results)

	save_results_df(
		out_dir, results_df, platform, sample_name, X_common_proteins, Y_rare_proteins
	)

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
	max_solutions,
    unique_reads_first_col_pos,
    unique_proteins_first_col_pos,
    x_and_y_proteins_denote_fractions::Bool = false,
)
    # Build a flat job list so we can thread over it safely
    jobs = Tuple[]
    
	for (platform, sample_name, expression_file, max_solution, unique_reads_file, unique_proteins_file) in zip(
		platforms, samples, expression_files, max_solutions, unique_reads_files, unique_proteins_files
	)
		push!(
			jobs, 
			(
				platform,
				sample_name,
				expression_file,
				max_solution,
				unique_reads_file,
				unique_proteins_file,
			)
		)
	end

    out = Vector{Union{Nothing, DataFrame}}(undef, length(jobs))

    Threads.@threads for idx in eachindex(jobs)
        platform,
        sample_name,
        expression_file,
        max_solution,
        unique_reads_file,
        unique_proteins_file = jobs[idx]

        # (Optional) avoid println here; it will interleave across threads.

        results_df = try
            process_one_sample(
                platform,
                sample_name,
                expression_file,
				max_solution,
                unique_proteins_file, 
				unique_proteins_first_col_pos,
                unique_reads_file, 
				unique_reads_first_col_pos,
                X_common_proteins, 
				Y_rare_proteins,
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

        out[idx] = results_df
    end

    # results_dfs = (df for df in out if df !== nothing)
    # results_df = isempty(results_dfs) ? DataFrame() : vcat(results_dfs...)

    # # If everything was skipped, return empty frames with a minimal schema
    # if nrow(results_df) == 0
    #     return DataFrame(), results_df
    # end

    # stats_df = combine(
    #     groupby(
    #         results_df,
    #         [:Platform, :IsSoftComparison, :XYProteinsDenoteFractions, :XCommonProteins, :YRareProteins, :Sample],
    #     ),
    #     :ActualXCommonProteins => first => :ActualXCommonProteins,
    #     :ActualYRareProteins => first => :ActualYRareProteins,
    #     :EditingSites => unique => :EditingSites,
    #     :EditingPercents => unique => :EditingPercents,
    #     :IsChimeric => sum => :NumOfChimericReads,
    #     :TotalCommonUniqueReads => first => :TotalCommonUniqueReads,
    #     :TotalRareUniqueReads => first => :TotalRareUniqueReads,
    #     :TotalCommonReads => first => :TotalCommonReads,
    #     :TotalRareReads => first => :TotalRareReads,
    # )
    # insertcols!(
    #     stats_df,
    #     "NumOfChimericReads",
    #     "%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.ActualYRareProteins,
    #     after = true,
    # )

    # return stats_df, results_df
end


function find_max_solution_file(expression_file)
	dir = dirname(expression_file)
	prefix = basename(expression_file) * ".MaxSolution_" # e.g. "GRIA.DistinctUniqueProteins.ExpressionLevels.EntropyConsidered.csv.MaxSolution_79"

	matches = filter(p -> isfile(p) && startswith(basename(p), prefix),
					readdir(dir; join=true))

	length(matches) == 1 || error("Expected exactly one matching file, found $(length(matches)) for expression file $expression_file with prefix $prefix: $(matches)")

	max_solution_file = only(matches)
end


function find_max_solution(expression_file)
	max_solution_file = find_max_solution_file(expression_file)
	# max_solution_str = max_solution_file[end-1:end]
	max_solution_str = rsplit(basename(max_solution_file), "_")[end]
	max_solution = parse(Int, max_solution_str)
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

max_solutions = find_max_solution.(expression_files)

unique_reads_first_col_pos = 9
unique_proteins_first_col_pos = 15

# X_common_proteins = [0.01, 0.03, 0.05]
X_common_proteins = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
Y_rare_proteins = [0.1, 0.2, 0.3]

x_and_y_proteins_denote_fractions = true

out_dir = "D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2"
# out_dir = "D.pealeii/MpileupAndTranscripts/JointChimericReadsAnalysis2Fixed"







# # TODO uncomment to test one sample
# i = 1
# # i = 5
# platform = platforms[i]
# sample_name = samples[i]
# expression_file = expression_files[i]
# max_solution = max_solutions[i]
# unique_reads_file = unique_reads_files[i]
# unique_proteins_file = unique_proteins_files[i]


# results_df = process_one_sample(
# 	platform,
# 	sample_name,
#     expression_file,
# 	max_solution,
#     unique_proteins_file, 
# 	unique_proteins_first_col_pos,
#     unique_reads_file, 
# 	unique_reads_first_col_pos,
# 	X_common_proteins, 
# 	Y_rare_proteins
# )




# stats_df, results_df = per_platform_stats_df(
#     X_common_proteins,
#     Y_rare_proteins,
#     platforms,
#     samples,
#     unique_reads_files,
#     unique_proteins_files,
#     expression_files,
# 	max_solutions,
#     unique_reads_first_col_pos,
#     unique_proteins_first_col_pos,
#     x_and_y_proteins_denote_fractions::Bool = false,
# )


per_platform_stats_df(
    X_common_proteins,
    Y_rare_proteins,
    platforms,
    samples,
    unique_reads_files,
    unique_proteins_files,
    expression_files,
	max_solutions,
    unique_reads_first_col_pos,
    unique_proteins_first_col_pos,
    x_and_y_proteins_denote_fractions,
)


if abspath(PROGRAM_FILE) == @__FILE__
	# don't run the plotting code when this file is included or run interactively, only when it's run directly
	exit(0)
end











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

