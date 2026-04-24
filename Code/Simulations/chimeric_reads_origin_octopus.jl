using DataFrames
using CSV
using Random # for MersenneTwister
using StatsBase  # for StatsBase.sample
using CairoMakie
using Format
using Transducers
using Base.Threads
using Dates


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
	results_df.ChimericUniqueReadsCombinations = map(
		v -> join(string.(first.(v), "-", last.(v)), ";"),
		results_df.ChimericUniqueReadsCombinations
	)

	mkpath(out_dir)  # make absolutely sure it exists (good under threads)
	
	soft_comparison_interfix_str = soft_comparison ? "SoftComparison" : "StrictComparison"
	out_file = out_dir * "/" * basename(replace(original_in_file, ".reads.snps.csv.gz" => ".chimeric_reads.$soft_comparison_interfix_str.csv.gz"))
	
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



function process_one_sample(
	snps_reads_file,
	out_dir,
	soft_comparison::Bool = false,
)

	reads_df = DataFrame(CSV.File(snps_reads_file, delim = "\t", types = Dict("Read" => String)))

	# isempty(reads_df) && return DataFrame() # return empty DataFrame if there are no reads
	if isempty(reads_df)
        # Still write an output so every input has both {Soft,Strict} files.
        empty_results_df = DataFrame(
            Platform = String[],
            Chrom = String[],
            UniqueRead = String[],
            IsChimeric = Bool[],
            ChimericUniqueReadsCombinations = String[],
            NumOfChimericCombinations = Int[],
            IsSoftComparison = Bool[],
        )
        save_results_df(empty_results_df, snps_reads_file, out_dir, soft_comparison)
        return empty_results_df
    end

	rename!(reads_df, "Sample" => "Chrom")

	gdf = groupby(reads_df, Not("Platform", "Chrom", "Read"))

	reads_df = combine(
		gdf,
		"Platform" => unique => "Platform",
		"Chrom" => unique => "Chrom",
		"Read" => (x -> join(x, ",")) => "Reads",
	)

	reads_df = select(
		reads_df,
		"Platform",
		"Chrom",
		"Reads",
		Not("Platform", "Chrom", "Reads"),
	)

	platform = reads_df.Platform[1]
	chrom = reads_df.Chrom[1]

	# TODO 19.4 chat change
	# results = []
	results = Vector{NamedTuple}(undef, 0)

	for x in 1:nrow(reads_df)

		read_x_row = reads_df[x, :]
		read_x = read_x_row[:Reads]
		read_x_snps_status_array = Array(read_x_row[Not("Platform", "Chrom", "Reads")])

		# read_x_snps_status_array'

		reads_other_than_x_df = reads_df[Not(x), :]
		reads_other_than_x_snps_status_array = Array(reads_other_than_x_df[:, Not("Platform", "Chrom", "Reads")])

		if soft_comparison
			M = softcomparison.(eachrow(reads_other_than_x_snps_status_array), Ref(read_x_snps_status_array))
			M = hcat(M...)'
		else
			M = reads_other_than_x_snps_status_array .== read_x_snps_status_array'
		end

		is_chimeric = are_there_chimeric_reads(M)
		all_chimeric_reads_indices = find_all_chimeric_reads(M)

		chimeric_unique_reads_combinations = [
			(reads_other_than_x_df[i, "Reads"], reads_other_than_x_df[j, "Reads"])
			for (i, j) in all_chimeric_reads_indices
		]

		result = (
			Platform = platform,
			Chrom = chrom,
			UniqueRead = read_x,
			IsChimeric = is_chimeric,
			ChimericUniqueReadsCombinations = chimeric_unique_reads_combinations,
			NumOfChimericCombinations = length(chimeric_unique_reads_combinations),
			IsSoftComparison = soft_comparison,
		)
		push!(results, result)

	end

	results_df = DataFrame(results)

	save_results_df(
		results_df,
		snps_reads_file,
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
	snps_reads_files,
	out_dir,
)	
	# TODO 19.4 chat change
	mkpath(out_dir)

	Threads.@threads for snps_reads_file in snps_reads_files
		for soft_comparison in [true, false]
			# TODO 19.4 chat change
			# process_one_sample(
			# 	snps_reads_file,
			# 	out_dir,
			# 	soft_comparison,
			# )
			try
                process_one_sample(snps_reads_file, out_dir, soft_comparison)
            catch err
                # log and continue
                mode = soft_comparison ? "SoftComparison" : "StrictComparison"
                log_file = joinpath(out_dir, basename(snps_reads_file) * "." * mode * ".error.log")
                open(log_file, "w") do io
                    println(io, "time\t", Dates.now())
                    println(io, "thread\t", threadid())
                    println(io, "input\t", snps_reads_file)
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
out_dir = base_dir * "/ChimericReadsFiles"

# IMPORTANT: don't wipe outputs while debugging.
# try
# 	rm(out_dir, recursive=true)
# catch IOError
# 	# do nothing if the directory does not exist
# end
# mkdir(out_dir)

mkpath(out_dir)

snps_reads_files = filter(
	f -> endswith(f, ".reads.snps.csv.gz"),
	readdir(reads_dir; join = true),
)

# For debugging, start small:
# snps_reads_files = snps_reads_files[1:10]


# process_samples(
# 	snps_reads_files,
# 	out_dir,
# )







# corrupt_chroms = ["comp178222_c0_seq2"]

# just a few of them, not all 223 of them
missing_chroms = [
	"comp143597_c0_seq1",
	"comp143722_c0_seq1",
	# "comp143834_c0_seq2",
	# "comp143840_c0_seq4",
	# "comp143866_c0_seq1",
]

all_chroms = [
	split(split(snps_reads_file, "/")[end], ".")[1]
	for snps_reads_file in snps_reads_files
]

good_chroms = setdiff(all_chroms, union(corrupt_chroms, missing_chroms))


test_chroms = vcat(
	good_chroms[1:2],
	corrupt_chroms,
	missing_chroms[1:2]
) 

snps_reads_files = [
	snps_reads_file
	for snps_reads_file in snps_reads_files
		for chrom in test_chroms
			if occursin(chrom, snps_reads_file)
]


process_samples(
	snps_reads_files,
	out_dir,
)







missing_chroms = """comp143597_c0_seq1
comp143722_c0_seq1
comp143834_c0_seq2
comp143840_c0_seq4
comp143866_c0_seq1
comp143867_c0_seq1
comp143875_c0_seq1
comp143905_c0_seq1
comp143934_c0_seq1
comp143937_c0_seq2
comp143943_c0_seq1
comp143965_c0_seq1
comp143979_c0_seq1
comp143989_c1_seq1
comp144014_c0_seq1
comp144057_c0_seq2
comp144109_c0_seq1
comp144117_c0_seq2
comp144126_c0_seq1
comp144130_c1_seq3
comp144143_c0_seq1
comp144147_c1_seq1
comp144164_c0_seq1
comp159552_c0_seq4
comp159592_c0_seq1
comp159682_c0_seq1
comp159694_c0_seq1
comp159699_c0_seq1
comp159711_c0_seq1
comp159741_c0_seq1
comp159747_c0_seq1
comp159761_c0_seq1
comp159763_c0_seq2
comp159776_c0_seq1
comp159779_c0_seq1
comp159828_c0_seq1
comp159843_c0_seq1
comp160218_c0_seq1
comp160228_c0_seq2
comp160254_c0_seq1
comp160303_c1_seq1
comp160326_c0_seq1
comp160334_c0_seq1
comp160335_c0_seq1
comp160427_c0_seq2
comp160444_c0_seq2
comp160566_c0_seq2
comp160588_c0_seq2
comp160717_c1_seq1
comp160723_c0_seq1
comp160731_c0_seq1
comp160740_c0_seq1
comp160745_c0_seq1
comp160761_c0_seq1
comp160774_c0_seq2
comp160785_c0_seq1
comp160826_c0_seq1
comp160841_c0_seq1
comp160858_c0_seq2
comp160871_c0_seq1
comp160872_c0_seq1
comp160939_c0_seq1
comp160990_c0_seq1
comp160996_c0_seq3
comp160998_c0_seq1
comp161022_c1_seq2
comp161040_c0_seq1
comp161056_c0_seq1
comp161108_c0_seq1
comp161110_c0_seq3
comp161129_c0_seq1
comp161170_c0_seq1
comp161172_c0_seq2
comp161283_c0_seq2
comp161286_c0_seq1
comp161297_c0_seq1
comp161348_c0_seq2
comp161356_c0_seq1
comp161371_c0_seq1
comp161400_c1_seq1
comp161408_c0_seq1
comp161439_c0_seq1
comp161442_c0_seq1
comp161452_c0_seq1
comp161507_c0_seq1
comp178223_c1_seq1
comp178238_c0_seq2
comp178239_c3_seq1
comp178244_c0_seq1
comp178244_c2_seq1
comp178244_c4_seq1
comp178246_c0_seq4
comp178248_c0_seq1
comp178264_c0_seq1
comp178266_c0_seq26
comp178267_c0_seq1
comp178274_c1_seq1
comp178290_c0_seq3
comp178291_c0_seq1
comp178292_c0_seq1
comp178306_c1_seq1
comp178310_c0_seq1
comp178322_c0_seq1
comp178324_c0_seq2
comp178329_c0_seq2
comp178331_c0_seq8
comp178334_c0_seq1
comp178355_c0_seq1
comp178370_c0_seq1
comp178381_c0_seq1
comp178389_c0_seq2
comp178413_c0_seq1
comp178419_c0_seq1
comp178419_c1_seq1
comp178421_c0_seq3
comp178427_c1_seq13
comp178450_c0_seq3
comp178451_c0_seq1
comp178457_c0_seq2
comp178478_c0_seq6
comp178483_c0_seq2
comp178485_c0_seq2
comp178489_c0_seq1
comp179461_c1_seq5
comp179463_c1_seq4
comp179463_c2_seq2
comp179466_c8_seq1
comp179469_c0_seq6
comp179470_c3_seq1
comp179475_c0_seq5
comp179478_c0_seq2
comp179480_c0_seq1
comp179480_c1_seq1
comp179486_c0_seq1
comp179489_c1_seq1
comp179490_c0_seq11
comp179495_c2_seq1
comp179497_c2_seq1
comp179498_c0_seq1
comp179510_c3_seq4
comp179510_c5_seq8
comp179511_c2_seq1
comp179515_c4_seq1
comp179520_c1_seq1
comp179527_c1_seq2
comp179531_c2_seq1
comp179569_c0_seq2
comp179572_c0_seq1
comp179575_c1_seq2
comp179579_c1_seq1
comp179581_c0_seq1
comp179582_c0_seq2
comp179591_c0_seq6
comp179595_c3_seq1
comp179596_c4_seq1
comp179600_c0_seq1
comp179604_c0_seq1
comp179615_c0_seq2
comp181659_c1_seq4
comp181663_c0_seq1
comp181664_c0_seq1
comp181664_c1_seq3
comp181665_c3_seq1
comp181669_c0_seq2
comp181672_c0_seq2
comp181673_c0_seq1
comp181677_c0_seq5
comp181679_c0_seq1
comp181686_c0_seq2
comp181687_c0_seq7
comp181689_c0_seq1
comp181691_c0_seq8
comp181694_c0_seq1
comp181697_c0_seq2
comp181698_c0_seq2
comp181700_c0_seq4
comp181704_c0_seq3
comp181711_c1_seq3
comp181711_c6_seq1
comp181711_c7_seq1
comp181713_c11_seq1
comp181713_c4_seq1
comp181713_c9_seq1
comp181718_c0_seq1
comp181719_c0_seq1
comp181723_c2_seq2
comp181723_c3_seq1
comp181728_c1_seq1
comp181734_c0_seq1
comp181734_c2_seq8
comp181739_c2_seq14
comp181744_c0_seq1
comp181748_c0_seq4
comp181754_c0_seq4
comp181757_c0_seq23
comp181764_c0_seq1
comp181765_c0_seq2
comp181766_c0_seq2
comp181769_c0_seq1
comp181775_c3_seq2
comp181776_c0_seq20
comp181792_c0_seq2
comp181810_c0_seq1
comp181820_c2_seq1
comp181822_c1_seq1
comp1818282_c0_seq1
comp181832_c0_seq1
comp181833_c0_seq2
comp181838_c0_seq1
comp181838_c1_seq1
comp181844_c1_seq1
comp181844_c7_seq5
comp72309_c0_seq1
comp72317_c1_seq1
comp72326_c0_seq1
comp72330_c0_seq1
comp72344_c0_seq1
comp72348_c0_seq1
comp72358_c0_seq1
comp72368_c0_seq1
comp72374_c0_seq1
comp72375_c0_seq1
comp72376_c0_seq1"""
missing_chroms = split(missing_chroms, "\n")






# 



# snps_reads_file = snps_reads_files[3]
# soft_comparison = true
# results_df = process_one_sample(snps_reads_file, soft_comparison)

# stats_df = per_platform_stats_df(snps_reads_files)
# describe(stats_df[stats_df.IsSoftComparison, "%OfChimericReads"])
# describe(stats_df[stats_df.IsSoftComparison .== false, "%OfChimericReads"])









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



# function plot_stats_df(
# 	X_common_proteins,
# 	Y_rare_proteins,
# 	stats_df,
# 	samples,
# 	yscale = identity,
# 	colors = Makie.wong_colors(),
# )
# 	samples_cat = Dict(sample => i for (i, sample) in enumerate(samples))
# 	xticks = (
# 		1:length(samples),
# 		samples,
# 	)

# 	for soft_comparison in [true, false]

# 		fig = Figure(size = (600, 600))
# 		axes = []
# 		for (i, x_common_proteins) in enumerate(X_common_proteins)
# 			for (j, y_rare_proteins) in enumerate(Y_rare_proteins)
# 				title = "Common: $(format(x_common_proteins, commas=true))\nRare: $(format(y_rare_proteins, commas=true))"
# 				ax = Axis(
# 					fig[i, j],
# 					xticks = xticks,
# 					xticklabelrotation = π / 4,
# 					subtitle = title,
# 					yscale = yscale,
# 				)
# 				push!(axes, ax)
# 				subdf = stats_df[stats_df.XCommonProteins.==x_common_proteins.&&stats_df.YRareProteins.==y_rare_proteins.&&stats_df.IsSoftComparison.==soft_comparison, :]
# 				xs = [samples_cat[sample] for sample in subdf.Sample]
# 				ys = subdf[!, "%OfChimericReads"]

# 				barplot!(xs, ys,
# 					color = colors[xs],
# 					# strokecolor = :black, 
# 					# strokewidth = 1
# 				)
# 			end
# 		end
# 		linkxaxes!(axes...)
# 		linkyaxes!(axes...)
# 		Label(fig[begin:end, 0], "% of chimeric reads", rotation = pi / 2)  # y axis title

# 		if soft_comparison
# 			title = "Soft comparison"
# 		else
# 			title = "Strict comparison"
# 		end
# 		Label(fig[0, begin+1:end], title, fontsize = 18)  # main title

# 		display(fig)
# 	end

# end




# plot_stats_df(
# 	X_common_proteins,
# 	Y_rare_proteins,
# 	pacbio_stats_df,
# 	samples,
# 	# log2
# )


# pacbio_stats_df[pacbio_stats_df.Sample .== "ADAR1" .&& pacbio_stats_df.IsSoftComparison, :]
# pacbio_stats_df[pacbio_stats_df.Sample .== "ADAR1" .&& pacbio_stats_df.IsSoftComparison .== false, :]

# pacbio_stats_df[pacbio_stats_df.Sample .== "IQEC1" .&& pacbio_stats_df.IsSoftComparison, :]
# pacbio_stats_df[pacbio_stats_df.Sample .== "IQEC1" .&& pacbio_stats_df.IsSoftComparison .== false, :]






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