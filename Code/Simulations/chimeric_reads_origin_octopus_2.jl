using DataFrames
using CSV
using Random # for MersenneTwister
using StatsBase  # for StatsBase.sample
using CairoMakie
using Format
using Transducers
using Base.Threads
using Dates
using Transducers


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


function find_all_chimeric_reads(B, E, debug::Bool = false)
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
    # As i increases (sorted_B nondecreasing), the set {E[j] <= B[i]} only grows,
    # so we can advance j monotonically.
	pairs = Vector{Tuple{Int, Int}}()
	j = 1
	Eⱼ = sorted_E[j]
	N = length(sorted_B) # same as length(sorted_E)

	@inbounds for i ∈ 1:N
		
		debug && println("i = $i, j = $j")
		
		Bᵢ = sorted_B[i]
		while j <= N && Eⱼ <= Bᵢ
			j += 1
			Eⱼ = sorted_E[j]
		end
		
		debug && println("j = $j")
		
		# All E indices in 1:(j-1) satisfy E <= B[i]
		
		debug && println("k ∈ 1:j-1 = $(collect(1:j-1))")
		
		Bᵢ_orig = idxs_B[i]
		for k ∈ 1:j-1
			push!(pairs, (Bᵢ_orig, idxs_E[k]))
		end
	end

	@assert length(unique(pairs)) == length(pairs) "All pairs should be unique, i.e. we should not have duplicate pairs of (r_common, r_common2)"

	filter!(
		pair -> pair[1] != pair[2], # filter out pairs of the same read (r_common, r_common)
		pairs
	)

	@assert Set(length.(Set.(pairs))) == Set(2) "All pairs should be of length 2, i.e. tuples of (r_common, r_common2)"

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


# extract_chimerizing_pairs(reads_other_than_x_df::DataFrame, all_chimeric_reads_indices) = [
# 	(reads_other_than_x_df[i, "UniqueRead"], reads_other_than_x_df[j, "UniqueRead"])
# 	for (i, j) ∈ all_chimeric_reads_indices
# ]

extract_chimerizing_pairs(reads_other_than_x_names::Vector{String}, all_chimeric_reads_indices) = [
	(reads_other_than_x_names[i], reads_other_than_x_names[j])
	for (i, j) ∈ all_chimeric_reads_indices
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



function validate_chimerizing_pairs_are_made_of_different_reads(unique_reads_file, chimerizing_pairs, soft_comparison)
	for (i, (r1, r2)) in enumerate(chimerizing_pairs)
		r1 == r2 && error(
			"Found a chimeric pair made of identical reads during soft_comparison = $soft_comparison, which should not happen. 
			Please check the input file: $unique_reads_file. The first such pair is: ($r1, $r2) (#$i pair)"
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


function save_results_df(
	results_df,
	original_in_file,
	out_dir,
	soft_comparison
)
	results_df = deepcopy(results_df)
	# results_df.ChimericUniqueReadsCombinations = map(
	# 	v -> join(string.(first.(v), "-", last.(v)), ";"),
	# 	results_df.ChimericUniqueReadsCombinations
	# )
	results_df.ChimericUniqueReadsCombinations = map(
		v -> join(string.(first.(v), ",", last.(v)), ";"),
		results_df.ChimericUniqueReadsCombinations
	)
	results_df.ChimerizingSitesIntersections = map(
		v -> join(string.(first.(v), ",", last.(v)), ";"),
		results_df.ChimerizingSitesIntersections
	)

	mkpath(out_dir)  # make absolutely sure it exists (good under threads)
	
	# soft_comparison_interfix_str = soft_comparison ? "SoftComparison" : "StrictComparison"
	soft_comparison_interfix_str = soft_comparison ? "soft_comparison" : "strict_comparison"
	# out_file = out_dir * "/" * basename(replace(original_in_file, ".reads.snps.csv.gz" => ".chimeric_reads.$soft_comparison_interfix_str.csv.gz"))
	out_file = out_dir * "/" * basename(replace(original_in_file, ".unique_reads.csv.gz" => ".chimeric_reads.$soft_comparison_interfix_str.csv.gz"))
	
	# Atomic write: write to temp then move into place.
    tmp_file = out_file * ".tmp." * string(getpid()) * "." * string(threadid())
	
	# # bump buffer to handle very long rows (large joined string)
	# CSV.write(out_file, results_df, delim = "\t"; compress = true, bufsize = 64 * 1024 * 1024)

	# bump buffer to handle very long rows (large joined string)
	# Write to temp first, then move into place (atomic within same filesystem).
    CSV.write(tmp_file, results_df, delim = "\t"; compress = true, bufsize = 64 * 1024 * 1024)
	mv(tmp_file, out_file; force = true)
    
	return out_file
end



function prepare_unique_reads_df(chrom, unique_reads_file, unique_reads_first_col_pos)
	unique_reads_df = DataFrame(
		CSV.File(
			unique_reads_file, 
			delim = "\t", 
			types = Dict(
				"UniqueRead" => String,
				"Reads" => String,
			)
		)
	)
	# rename!(unique_reads_df, "Transcript" => "Gene")
	
	#= #todo by retaining "Samples" and corresponding "Reads" columns, 
	we can later more correctly only compare reads that are present in the same sample, 
	instead of comparing all reads across all samples
	(which is inflating the number of chimeric reads because chimeras are only created in the PCR per sample) =#
	
	# select!(unique_reads_df, vcat(["Gene", "UniqueRead", "NumOfReads"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	# insertcols!(unique_reads_df, 4, "Chrom" => chrom)
	# new_unique_reads_first_col_pos = 4

	select!(unique_reads_df, vcat(["UniqueRead", "NumOfReads"], names(unique_reads_df)[unique_reads_first_col_pos:end]))
	insertcols!(unique_reads_df, 1, "Chrom" => chrom)
	new_unique_reads_first_col_pos = 4
	return unique_reads_df, new_unique_reads_first_col_pos
end



function process_one_sample(
	platform,
	chrom,
	unique_reads_file,
	unique_reads_first_col_pos,
	out_dir,
	soft_comparison::Bool = false,
)

	unique_reads_df, new_unique_reads_first_col_pos = prepare_unique_reads_df(
		chrom, unique_reads_file, unique_reads_first_col_pos
	)

	unique_reads_statuses_df = unique_reads_df[:, new_unique_reads_first_col_pos:end]
	nrow(unique(unique_reads_statuses_df)) == nrow(unique_reads_statuses_df) || error(
		"There are duplicate rows in the unique reads editing statuses, which should not happen. Please check the input file: $unique_reads_file"
	)
	editing_sites = parse.(Int, names(unique_reads_statuses_df))
	editing_percents = [
		100 * count(==(1), col) / count(!=(-1), col)
		for col ∈ eachcol(unique_reads_statuses_df)
	]


	# isempty(reads_df) && return DataFrame() # return empty DataFrame if there are no reads
	if isempty(unique_reads_df)
        # Still write an output so every input has both {Soft,Strict} files.
        empty_results_df = DataFrame(
            Platform = String[],
            Chrom = String[],
			EditingSites = Int[],
			EditingPercents = Float64[],
			IsSoftComparison = Bool[],
            UniqueRead = String[],
			NumOfReads = Int[],
            IsChimeric = Bool[],
            NumOfChimericCombinations = Int[],
			ChimericUniqueReadsCombinations = Tuple{String,String}[],
            ChimerizingSitesIntersections = Tuple{Int,Int}[],
        )
        save_results_df(empty_results_df, unique_reads_file, out_dir, soft_comparison)
        return empty_results_df
    end

	results = Vector{NamedTuple}(undef, 0)

	for x ∈ 1:nrow(unique_reads_df)

		# x = 1 # todo comment-out, this is for testing

		read_x_row = unique_reads_df[x, :]
		read_x = read_x_row[:UniqueRead]
		read_x_editing_status_array = Array(read_x_row[Not("Chrom", "UniqueRead", "NumOfReads")])
		read_x_num_of_reads = read_x_row[:NumOfReads]

		# read_x_snps_status_array'

		reads_other_than_x_df = unique_reads_df[Not(x), :]
		reads_other_than_x_editing_status_array = Array(reads_other_than_x_df[:, Not("Chrom", "UniqueRead", "NumOfReads")])

		if soft_comparison
			M = softcomparison.(eachrow(reads_other_than_x_editing_status_array), Ref(read_x_editing_status_array))
			M = hcat(M...)'
		else
			M = reads_other_than_x_editing_status_array .== read_x_editing_status_array'
		end

		# is_chimeric = are_there_chimeric_reads(M)
		# all_chimeric_reads_indices = find_all_chimeric_reads(M)

		B, E = find_B_and_E_for_M(M)
		is_chimeric = are_there_chimeric_reads(B, E)
		all_chimeric_reads_indices = find_all_chimeric_reads(B, E)
		chimerizing_pairs = extract_chimerizing_pairs(reads_other_than_x_df, all_chimeric_reads_indices)
		chimeric_reads_sites_intersection = extract_chimeric_reads_sites_intersection(all_chimeric_reads_indices, B, E)

		# validate_chimerizing_pairs_are_made_of_different_reads(unique_reads_file, soft_comparison, chimerizing_pairs)

		# i = 152
		# all_chimeric_reads_indices[i]
		# chimerizing_pairs[i]
		# chimeric_reads_sites_intersection[i]
		# B[8]
		# E[1]

		result = (
			Platform = platform,
			Chrom = chrom,
			EditingSites = editing_sites,
			EditingPercents = editing_percents,
			IsSoftComparison = soft_comparison,
			UniqueRead = read_x,
			NumOfReads = read_x_num_of_reads,
			IsChimeric = is_chimeric,
			NumOfChimericCombinations = length(chimerizing_pairs),
			ChimericUniqueReadsCombinations = chimerizing_pairs,
			ChimerizingSitesIntersections = chimeric_reads_sites_intersection,
		)
		push!(results, result)

	end

	results_df = DataFrame(results)

	save_results_df(
		results_df,
		unique_reads_file,
		out_dir,
		soft_comparison
	)

	return results_df
end



function process_one_sample_strict_comparison(
	platform,
	chrom,
	unique_reads_file,
	unique_reads_first_col_pos,
	out_dir,
)
	soft_comparison = false

	unique_reads_df, new_unique_reads_first_col_pos = prepare_unique_reads_df(
		chrom, unique_reads_file, unique_reads_first_col_pos
	)

	unique_reads_statuses_df = unique_reads_df[:, new_unique_reads_first_col_pos:end]
	nrow(unique(unique_reads_statuses_df)) == nrow(unique_reads_statuses_df) || error(
		"There are duplicate rows in the unique reads editing statuses, which should not happen. Please check the input file: $unique_reads_file"
	)
	editing_sites_names = parse.(Int, names(unique_reads_statuses_df))
	num_of_editing_sites = length(editing_sites_names)
	editing_percents = [
		100 * count(==(1), col) / count(!=(-1), col)
		for col ∈ eachcol(unique_reads_statuses_df)
	]

	# isempty(reads_df) && return DataFrame() # return empty DataFrame if there are no reads
	if isempty(unique_reads_df)
        # Still write an output so every input has both {Soft,Strict} files.
        empty_results_df = DataFrame(
            Platform = String[],
            Chrom = String[],
			NumOfEditingSites = Int[],
			EditingSites = Int[],
			EditingPercents = Float64[],
			IsSoftComparison = Bool[],
            UniqueRead = String[],
			NumOfReads = Int[],
            IsChimeric = Bool[],
            NumOfChimericCombinations = Int[],
			ChimericUniqueReadsCombinations = Tuple{String,String}[],
            ChimerizingSitesIntersections = Tuple{Int,Int}[],
        )
        save_results_df(empty_results_df, unique_reads_file, out_dir, soft_comparison)
        return empty_results_df
    end

	results = Vector{NamedTuple}(undef, 0)

	for x ∈ 1:nrow(unique_reads_df)

		# x = 1

		read_x_row = unique_reads_df[x, :]
		read_x_name = read_x_row[:UniqueRead]
		read_x_editing_status_array = Array(read_x_row[Not("Chrom", "UniqueRead", "NumOfReads")])
		read_x_num_of_reads = read_x_row[:NumOfReads]

		# read_x_snps_status_array'

		reads_other_than_x_df = unique_reads_df[Not(x), :]
		reads_other_than_x_names = reads_other_than_x_df[:, :UniqueRead]
		reads_other_than_x_editing_status_array = Array(reads_other_than_x_df[:, Not("Chrom", "UniqueRead", "NumOfReads")])

		# if soft_comparison
		# 	M = softcomparison.(eachrow(reads_other_than_x_editing_status_array), Ref(read_x_editing_status_array))
		# 	M = hcat(M...)'
		# else
		# 	M = reads_other_than_x_editing_status_array .== read_x_editing_status_array'
		# end
		M = reads_other_than_x_editing_status_array .== read_x_editing_status_array'

		@assert !any(sum.(eachrow(M)) .== size(M, 2)) "There is another unique read that is identical to the considered unique read $read_x_name, which should not happen. Please check the input file: $unique_reads_file. The first such read is: $(reads_other_than_x_df[findfirst(sum.(eachrow(M)) .== size(M, 2)), :UniqueRead])"

		B, E = find_B_and_E_for_M(M)
		# is_chimeric = are_there_chimeric_reads(B, E)
		all_chimeric_reads_indices = find_all_chimeric_reads(B, E)
		is_chimeric = length(all_chimeric_reads_indices) > 0
		# chimerizing_pairs = extract_chimerizing_pairs(reads_other_than_x_df, all_chimeric_reads_indices)
		chimerizing_pairs = extract_chimerizing_pairs(reads_other_than_x_names, all_chimeric_reads_indices)
		chimeric_reads_sites_intersection = extract_chimeric_reads_sites_intersection(all_chimeric_reads_indices, B, E)

		validate_chimerizing_pairs_are_made_of_different_reads(
			unique_reads_file, chimerizing_pairs, soft_comparison
		)

		result = (
			Platform = platform,
			Chrom = chrom,
			NumOfEditingSites = num_of_editing_sites,
			EditingSites = editing_sites_names,
			EditingPercents = editing_percents,
			IsSoftComparison = soft_comparison,
			UniqueRead = read_x_name,
			NumOfReads = read_x_num_of_reads,
			IsChimeric = is_chimeric,
			NumOfChimericCombinations = length(chimerizing_pairs),
			ChimericUniqueReadsCombinations = chimerizing_pairs,
			ChimerizingSitesIntersections = chimeric_reads_sites_intersection,
		)
		push!(results, result)

	end

	results_df = DataFrame(results)

	save_results_df(
		results_df,
		unique_reads_file,
		out_dir,
		soft_comparison
	)

	return results_df
end




function process_one_sample_soft_comparison(
	platform,
	chrom,
	unique_reads_file,
	unique_reads_first_col_pos,
	out_dir,
)
	soft_comparison = true

	unique_reads_df, new_unique_reads_first_col_pos = prepare_unique_reads_df(
		chrom, unique_reads_file, unique_reads_first_col_pos
	)

	unique_reads_statuses_df = unique_reads_df[:, new_unique_reads_first_col_pos:end]
	nrow(unique(unique_reads_statuses_df)) == nrow(unique_reads_statuses_df) || error(
		"There are duplicate rows in the unique reads editing statuses, which should not happen. Please check the input file: $unique_reads_file"
	)
	editing_sites_names = parse.(Int, names(unique_reads_statuses_df))
	num_of_editing_sites = length(editing_sites_names)
	editing_percents = [
		100 * count(==(1), col) / count(!=(-1), col)
		for col ∈ eachcol(unique_reads_statuses_df)
	]


	
	# # 1. Allow your DataFrame to hold 'missing' values
	# # explicitly_missing_unique_reads_df = copy(unique_reads_df)
	# # allowmissing!(explicitly_missing_unique_reads_df)
	# explicitly_missing_unique_reads_df = allowmissing(unique_reads_df)
	# # 2. Perform the replacement (notice 'replace' instead of 'replace!')
	# explicitly_missing_unique_reads_df = mapcols(c -> replace(c, -1 => missing), explicitly_missing_unique_reads_df)

	# gdf = groupby(
	# 	explicitly_missing_unique_reads_df,
	# 	names(unique_reads_statuses_df)[new_unique_reads_first_col_pos:end]
	# )
	# gdf[1]
	# gdf[2]


	# isempty(reads_df) && return DataFrame() # return empty DataFrame if there are no reads
	if isempty(unique_reads_df)
        # Still write an output so every input has both {Soft,Strict} files.
        empty_results_df = DataFrame(
            Platform = String[],
            Chrom = String[],
			NumOfEditingSites = Int[],
			EditingSites = Int[],
			EditingPercents = Float64[],
			IsSoftComparison = Bool[],
            UniqueRead = String[],
			NumOfReads = Int[],
            IsChimeric = Bool[],
            NumOfChimericCombinations = Int[],
			ChimericUniqueReadsCombinations = Tuple{String,String}[],
            ChimerizingSitesIntersections = Tuple{Int,Int}[],
        )
        save_results_df(empty_results_df, unique_reads_file, out_dir, soft_comparison)
        return empty_results_df
    end

	results = Vector{NamedTuple}(undef, 0)

	for x ∈ 1:nrow(unique_reads_df)

		# x = 1

		read_x_row = unique_reads_df[x, :]
		read_x_name = read_x_row[:UniqueRead]

		read_x_editing_status_array = Array(read_x_row[Not("Chrom", "UniqueRead", "NumOfReads")])

		# missing_sites_in_read_x_bools = map(
		# 	site -> site == -1,
		# 	read_x_row[Not("Chrom", "UniqueRead", "NumOfReads")]
		# );
		# missing_sites_in_read_x = [
		# 	site
		# 	for (site, is_missing) in zip(
		# 		names(read_x_row[Not("Chrom", "UniqueRead", "NumOfReads")]),
		# 		missing_sites_in_read_x_bools
		# 	) if is_missing
		# ]
		# cols_to_omit_to_obatin_editing_status_array_considering_read_x = vcat(
		# 	["Chrom", "UniqueRead", "NumOfReads"],
		# 	missing_sites_in_read_x
		# )

		# read_x_editing_status_array = Array(
		# 	read_x_row[Not(cols_to_omit_to_obatin_editing_status_array_considering_read_x)]
		# )

		read_x_num_of_reads = read_x_row[:NumOfReads]

		# read_x_snps_status_array'

		reads_other_than_x_df = unique_reads_df[Not(x), :]
		reads_other_than_x_names = reads_other_than_x_df[:, :UniqueRead]
		reads_other_than_x_editing_status_array = Array(reads_other_than_x_df[:, Not("Chrom", "UniqueRead", "NumOfReads")])
		# reads_other_than_x_editing_status_array = Array(
		# 	reads_other_than_x_df[:, Not(cols_to_omit_to_obatin_editing_status_array_considering_read_x)]
		# )

		# if soft_comparison
		# 	M = softcomparison.(eachrow(reads_other_than_x_editing_status_array), Ref(read_x_editing_status_array))
		# 	M = hcat(M...)'
		# else
		# 	M = reads_other_than_x_editing_status_array .== read_x_editing_status_array'
		# end
		M = softcomparison.(eachrow(reads_other_than_x_editing_status_array), Ref(read_x_editing_status_array))
		M = hcat(M...)'

		
		
		reads_softly_identical_to_compared = reads_other_than_x_names[(sum.(eachrow(M)) .== size(M, 2))]

	
		B, E = find_B_and_E_for_M(M)
		# is_chimeric = are_there_chimeric_reads(B, E)
		all_chimeric_reads_indices = find_all_chimeric_reads(B, E)
		# chimerizing_pairs = extract_chimerizing_pairs(reads_other_than_x_df, all_chimeric_reads_indices)
		chimerizing_pairs = extract_chimerizing_pairs(reads_other_than_x_names, all_chimeric_reads_indices)

		# todo check if any (r1, r2) pairs ∈ chimerizing_pairs also (r1, r2) ∈ reads_softly_identical_to_compared, and if so, we need to filter them out (and from all_chimeric_reads_indices accordingly)

		chimeric_reads_sites_intersection = extract_chimeric_reads_sites_intersection(all_chimeric_reads_indices, B, E)

		validate_chimerizing_pairs_are_made_of_different_reads(unique_reads_file, soft_comparison, chimerizing_pairs)

		

		i = 152
		all_chimeric_reads_indices[i]
		chimerizing_pairs[i]
		chimeric_reads_sites_intersection[i]
		B[8]
		E[1]

		result = (
			Platform = platform,
			Chrom = chrom,
			NumOfEditingSites = num_of_editing_sites,
			EditingSites = editing_sites_names,
			EditingPercents = editing_percents,
			IsSoftComparison = soft_comparison,
			UniqueRead = read_x_name,
			NumOfReads = read_x_num_of_reads,
			IsChimeric = is_chimeric,
			NumOfChimericCombinations = length(chimerizing_pairs),
			ChimericUniqueReadsCombinations = chimerizing_pairs,
			ChimerizingSitesIntersections = chimeric_reads_sites_intersection,
		)
		push!(results, result)

	end

	results_df = DataFrame(results)

	save_results_df(
		results_df,
		unique_reads_file,
		out_dir,
		soft_comparison
	)

	return results_df
end


# function per_platform_stats_df(
# 	snps_reads_files,
# )
# 	results_dfs = tcollect(
# 		process_one_sample(
# 			snps_reads_file,
# 			soft_comparison,
# 		)
# 		for snps_reads_file in snps_reads_files
# 		for soft_comparison in [true, false]
# 	)
# 	results_dfs = [df for df in results_dfs if !isempty(df)]

# 	results_df = vcat(results_dfs...)

# 	stats_df = combine(
# 		groupby(results_df, [:Platform, :IsSoftComparison, :Chrom]),
# 		:IsChimeric => sum => :NumOfChimericReads,
# 		nrow => :NumOfUniqueReads,
# 	)
# 	insertcols!(stats_df, "%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.NumOfUniqueReads)

# 	return stats_df
# end 


function process_samples(
	platform,
	chroms,
	unique_reads_files,
	unique_reads_first_col_pos,
	out_dir,
)
	mkpath(out_dir)

	# Threads.@threads for (chrom, unique_reads_file) ∈ zip(chroms, unique_reads_files)
	Threads.@threads for (chrom, unique_reads_file) ∈ collect(zip(chroms, unique_reads_files))
		for soft_comparison in [true, false]
			try
                process_one_sample(
					platform,
					chrom,
					unique_reads_file,
					unique_reads_first_col_pos,
					out_dir,
					soft_comparison,
				)
            catch err
                # log and continue
                mode = soft_comparison ? "SoftComparison" : "StrictComparison"
                log_file = joinpath(out_dir, basename(unique_reads_file) * "." * mode * ".error.log")
                open(log_file, "w") do io
                    println(io, "time\t", Dates.now())
                    println(io, "thread\t", threadid())
                    println(io, "input\t", unique_reads_file)
                    println(io, "mode\t", mode)
                    println(io, "error\t", sprint(showerror, err))
                    println(io, "backtrace:")
                    show(io, Base.catch_backtrace())
				end
            end
		end 
	end
end
# end




# ---- runner ----
base_dir = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.PooledSamples"
reads_dir = base_dir * "/ReadsFiles"
# out_dir = base_dir * "/ChimericReadsFiles"
out_dir = base_dir * "/ChimericReadsFiles2"

platform = "WholeTranscriptomePacBio"
unique_reads_first_col_pos = 10

# IMPORTANT: don't wipe outputs while debugging.
# try
# 	rm(out_dir, recursive=true)
# catch IOError
# 	# do nothing if the directory does not exist
# end
# mkdir(out_dir)

mkpath(out_dir)

unique_reads_files = filter(
	f -> endswith(f, ".unique_reads.csv.gz"),
	readdir(reads_dir; join = true),
)
chroms = [
	split(basename(f), ".")[1] for f in unique_reads_files
]

# # For debugging, start small:
# i = 10
# unique_reads_files = unique_reads_files[1:i]
# chroms = chroms[1:i]

unique_reads_file = unique_reads_files[1]
chrom = chroms[1]



process_samples(
	platform,
	chroms,
	unique_reads_files,
	unique_reads_first_col_pos,
	out_dir,
)




#  Parse back the serialized pair-lists:
#   "a,b;c,d"  -> [(a,b), (c,d)]
# Empty string / missing -> empty Vector{Tuple{String,String}}
parse_pairs_str(s) = begin
    if s === missing
        return Tuple{String,String}[]
    end
    str = String(s)
    isempty(str) && return Tuple{String,String}[]

    parts = split(str, ';'; keepempty=false)
    out = Tuple{String,String}[]

    for p in parts
        isempty(p) && continue
        ab = split(p, ','; limit=2)
        if length(ab) != 2 || isempty(ab[1]) || isempty(ab[2])
            # If you prefer hard-fail instead of skipping, replace with `error(...)`
            @warn "Bad pair token while parsing pairs" token=p raw=str
            continue
        end
        push!(out, (String(ab[1]), String(ab[2])))
    end

    return out
end

# If you want these as Tuple{Int,Int} instead of strings, use this:
parse_int_pairs_str(s) = begin
    if s === missing
        return Tuple{Int,Int}[]
    end
    str = String(s)
    isempty(str) && return Tuple{Int,Int}[]

    parts = split(str, ';'; keepempty=false)
    out = Tuple{Int,Int}[]

    for p in parts
        isempty(p) && continue
        ab = split(p, ','; limit=2)
        if length(ab) != 2 || isempty(ab[1]) || isempty(ab[2])
            @warn "Bad int-pair token while parsing pairs" token=p raw=str
            continue
        end
        try
            push!(out, (parse(Int, ab[1]), parse(Int, ab[2])))
        catch
            @warn "Non-integer token while parsing int pairs" token=p raw=str
        end
    end

    return out
end




function read_results_file(results_file)
	results_df = DataFrame(
	CSV.File(
			results_file, 
			delim = "\t", 
			types = Dict(
				"Platform" => String,
				"Chrom" => String,
				"UniqueRead" => String,
				"IsChimeric" => Bool,
				"NumOfChimericCombinations" => Int,
				"ChimericUniqueReadsCombinations" => String,
				"ChimerizingSitesIntersections" => String,
				"IsSoftComparison" => Bool,
			)
		)
	)
	try
		results_df.ChimericUniqueReadsCombinations = map(parse_pairs_str, results_df.ChimericUniqueReadsCombinations)
		results_df.ChimerizingSitesIntersections = map(parse_int_pairs_str, results_df.ChimerizingSitesIntersections)
	catch err
		# @warn "Error parsing pairs in results file" file=results_file error=sprint
		error("Failed parsing $results_file: $(sprint(showerror, err))")
	end
	
	return results_df
end



results_files = filter(
	f -> occursin(".chimeric_reads.", f),
	readdir(out_dir; join = true),
)



results_dfs = tcollect(
	read_results_file(results_file)
	for results_file in results_files
)

results_df = vcat(results_dfs...)

length(unique(results_df.Chrom))


@assert length(groupby(results_df, [:Chrom, :IsSoftComparison])) == length(chroms) * 2


stats_df = combine(
	groupby(
		results_df,
		[:Platform, :Chrom, :IsSoftComparison,],
	),
	# :EditingSites => unique => :EditingSites,
	# :EditingPercents => unique => :EditingPercents,
	:IsChimeric => sum => :NumOfChimericReads,
	nrow => :NumOfReads,
)
insertcols!(
	stats_df, 
	"NumOfChimericReads", 
	"%OfChimericReads" => 100 .* stats_df.NumOfChimericReads ./ stats_df.NumOfReads, 
	after=true
)



@assert sum(stats_df.NumOfReads) == nrow(results_df)


soft_comparisons = [false, true]


fig = Figure(size = (600, 500))
axes = Axis[]
row_1_axes = []
row_2_axes = []
xticks = 0:20:100
for (i, soft_comparison) in enumerate(soft_comparisons)
	sub_df = filter(row -> row.IsSoftComparison == soft_comparison, stats_df)
	# plot histograms of % chimeric reads per gene
	ax = Axis(
		fig[1, i], 
		title = soft_comparison ? "Soft comparison" : "Strict comparison", 
		xticks = xticks,
		ylabel = i == 1 ? "Genes" : "",
		yticks = 0:0.2:1.0
	)
	hist!(
		ax, sub_df[!, "%OfChimericReads"],
		normalization = :probability,
		# color = :blue
	)
	# ecdfplot!(
	# 	ax, sub_df[!, "%OfChimericReads"],
	# 	linecolor = :red,   # <- forces the ECDF line
	# 	color = :red
	# 	# normalization = :probability
	# )
	p = ecdfplot!(ax, sub_df[!, "%OfChimericReads"])
	# ecdfplot! may return a plot object or a vector/tuple of them; handle both:
	if p isa AbstractVector || p isa Tuple
		for q in p
			q.color = :red
		end
	else
		p.color = :red
	end
	push!(axes, ax)
	push!(row_1_axes, ax)
	# plot scatter of number of chimeric reads vs % chimeric reads per gene
	ax = Axis(
		fig[2, i], 
		xlabel = "Chimeric unique reads [%]", 
		xticks = xticks,
		ylabel = i == 1 ? "Unique reads" : "",
		yscale = log10
	)
	scatter!(
		ax, sub_df[!, "%OfChimericReads"], sub_df[!, "NumOfChimericReads"],
	)
	push!(axes, ax)
	push!(row_2_axes, ax)
end
linkyaxes!(row_1_axes...)
linkyaxes!(row_2_axes...)
display(fig)


