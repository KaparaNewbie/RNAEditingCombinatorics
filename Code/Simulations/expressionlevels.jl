using ArgParse
using CSV
using DataFrames
using StatsBase  # for StatsBase.sample
using BenchmarkTools
using BioSequences # for BioSequences.toAAset
import Base.Threads.@spawn
using Base.Threads
using ThreadsX # for ThreadsX.Set, 
using Transducers  # for tcollect
using SymmetricFormats # for SymmetricPacked which creates a symmetric matrix from a triangular array
using BioAlignments # for SubstitutionMatrix
using LinearAlgebra
using DataStructures # for Counter
using CairoMakie
using Logging, LoggingExtras


include(joinpath(@__DIR__, "consts.jl")) # for `AA_groups` (polar/non-polar/positive/negative)
include(joinpath(@__DIR__, "issimilar.jl")) # for issimilar
include(joinpath(@__DIR__, "timeformatters.jl"))


function prepare_distinctdf(
	distinctfile, delim, innerdelim, truestrings, falsestrings,
)
	@info "$(loggingtime())\tprepare_distinctdf" distinctfile
	# distinctdf = DataFrame(CSV.File(distinctfile; delim, truestrings, falsestrings))
	distinctdf = DataFrame(CSV.File(distinctfile; delim, truestrings, falsestrings, types = Dict("UniqueSamples" => String, "AvailableReads" => String)))
	# distinctdf[!, "UniqueSamples"] = InlineString.(distinctdf[!, "UniqueSamples"])

	transform!(distinctdf, :UniqueSamples => (x -> split.(x, innerdelim)) => :UniqueSamples)
	if "AvailableReads" ∈ names(distinctdf)
		transform!(distinctdf, :AvailableReads => (x -> split.(x, innerdelim)) => :AvailableReads)
	end
	distinctdf[!, "Index"] = collect(1:size(distinctdf, 1))
	return distinctdf
end


toAAset(x, innerdelim) = Set(map(aa -> convert(AminoAcid, only(aa)), split(x, innerdelim)))


function prepare_readsdf(readsfile, delim)
	@info "$(loggingtime())\tprepare_readsdf" readsfile delim
	readsdf = DataFrame(
		CSV.File(
			readsfile,
			delim = delim,
			types = Dict("Read" => String),
			select = ["Gene", "Read", "AmbigousPositions", "EditedPositions", "UneditedPositions"],
		),
	)
	rename!(readsdf, "AmbigousPositions" => "NAPositions")
	return readsdf
end


function prepare_allprotsdf!(
	allprotsfile, delim, innerdelim,
	truestrings, falsestrings, firstcolpos,
	readsdf,
)
	@info "$(loggingtime())\tprepare_allprotsdf!" allprotsfile

	df1 = DataFrame(CSV.File(allprotsfile, delim = delim, select = collect(1:firstcolpos-1), types = Dict("Protein" => String, "Reads" => String)))
	df1[!, "Protein"] = InlineString.(df1[!, :Protein])
	# make sure columns of AAs containing only Ts aren't parsed as boolean columns
	df2 = DataFrame(CSV.File(allprotsfile, delim = delim, drop = collect(1:firstcolpos-1), types = String))
	df2 = toAAset.(df2, innerdelim)

	allprotsdf = hcat(df1, df2)

	transform!(allprotsdf, :Reads => (x -> split.(x, innerdelim)) => :Reads)


	# calculate the mean na positions per each protein's reads - the current stats are inaccurate
	# readsdf = prepare_readsdf(readsfile, delim)
	mean_na_positions_per_prot = []
	for row in eachrow(allprotsdf)
		rows_reads = row["Reads"]
		na_positions_of_reads = readsdf[readsdf[!, "Read"].∈Ref(rows_reads), "NAPositions"]
		mean_na_positions_of_reads = mean(na_positions_of_reads)
		push!(mean_na_positions_per_prot, mean_na_positions_of_reads)
	end
	allprotsdf[!, "AmbigousPositions"] .= mean_na_positions_per_prot
	rename!(allprotsdf, "AmbigousPositions" => "MeanNAPositions")

	insertcols!(allprotsdf, firstcolpos, :Index => 1:size(allprotsdf, 1))
	firstcolpos += 1
	insertcols!(allprotsdf, firstcolpos, :AdditionalEqualSupportingReads => 0.0)
	firstcolpos += 1
	insertcols!(allprotsdf, firstcolpos, :AdditionalWeightedSupportingReads => 0.0)
	firstcolpos += 1

	insertcols!(
		allprotsdf,
		firstcolpos,
		:AdditionalEqualSupportingReadsContributionPerProtein => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalWeightedSupportingReadsContributionPerProtein => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalSupportingReadsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalSupportingProteinsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalSupportingProteinsDistances => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalSupportingProteinsMeanNAPositions => [[] for _ ∈ 1:size(allprotsdf, 1)],
	)
	firstcolpos += 6

	return allprotsdf, firstcolpos
end


const UInts = [UInt8, UInt16, UInt32, UInt64, UInt128]


"""Return the smallest unsigned int type that can hold `maximumval`."""
function smallestuint(maximumval)
	for uint ∈ UInts
		if maximumval <= typemax(uint)
			return uint
		end
		return UInts[end]
	end
end


"""Return pairs of sets of AAs with no two possible AAs that share a similar classification according to `aagroups`."""
function finddistinctAAsets(
	AAsets::Set{Set{AminoAcid}},
	aagroups::Dict{AminoAcid, String},
)
	@info "$(loggingtime())\tfinddistinctAAsets" AAsets aagroups
	ThreadsX.Set(
		[
		(x, y)
		for x ∈ AAsets for y ∈ AAsets
		if !anysimilarity(x, y, aagroups)
	]
	)
end


"""
Return pairs of sets of AAs with no two possible AAs such that their substitution score 
according to `substitutionmatrix` is `>`/`≥`/`<`/`≤` `similarityscorecutoff`. 
The specific comparison (e.g., `≥`) is determined by `similarityvalidator`.
"""
function finddistinctAAsets(
	AAsets::Set{Set{AminoAcid}},
	substitutionmatrix::SubstitutionMatrix{AminoAcid, Int64}, similarityscorecutoff::Int64, similarityvalidator::Function,
)
	@info "$(loggingtime())\tfinddistinctAAsets" AAsets substitutionmatrix similarityscorecutoff similarityvalidator
	ThreadsX.Set(
		[
		(x, y)
		for x ∈ AAsets for y ∈ AAsets
		if !anysimilarity(x, y, substitutionmatrix, similarityscorecutoff, similarityvalidator)
	]
	)
end


"""Return pairs of sets of AAs with no two possible shared AAs."""
function finddistinctAAsets(AAsets::Set{Set{AminoAcid}})
	@info "$(loggingtime())\tfinddistinctAAsets" AAsets
	ThreadsX.Set(
		[(x, y)
		 for x ∈ AAsets
		 for y ∈ AAsets
		 if x ∩ y == Set()]
	)
end


"""
	distances(M, aagroups)

Create a symmetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two contain no possible amino acids `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, 
such that both `AAᵦ` and `AAᵧ` share the same classification in `aagroups`.
"""
function distances(M::Matrix{Set{AminoAcid}}, aagroups::Dict{AminoAcid, String})
	@info "$(loggingtime())\tdistances" aagroups
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	distinctAAsets = finddistinctAAsets(AAsets, aagroups)
	return distances(M, distinctAAsets)
end


"""
	distances(M, substitutionmatrix, similarityscorecutoff, similarityvalidator)

Create a symmetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acis `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, 
such that their substitution score according to `substitutionmatrix` is `>`/`≥`/`<`/`≤` `similarityscorecutoff`. 
The specific comparison (e.g., `≥`) is determined by `similarityvalidator`.
"""
function distances(
	M::Matrix{Set{AminoAcid}},
	substitutionmatrix::SubstitutionMatrix{AminoAcid, Int64}, similarityscorecutoff::Int64, similarityvalidator::Function,
)
	@info "$(loggingtime())\tdistances" substitutionmatrix similarityscorecutoff similarityvalidator
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	distinctAAsets = finddistinctAAsets(AAsets, substitutionmatrix, similarityscorecutoff, similarityvalidator)
	return distances(M, distinctAAsets)
end


"""
	distances(M)

Create a symmetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acid.
"""
function distances(M::Matrix{Set{AminoAcid}})
	@info "$(loggingtime())\tdistances"
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	distinctAAsets = finddistinctAAsets(AAsets)
	return distances(M, distinctAAsets)
end


# M = [
# 	Set([AA_S]) Set([AA_S]) Set([AA_S])
# 	Set([AA_S]) Set([AA_S]) Set([AA_S])
# ]
# AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
# distinctAAsets = finddistinctAAsets(AAsets)

# distances(M, distinctAAsets)

"""
	distances(M, distinctAAsets)

Create a symmetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acid according to `distinctAAsets`.
"""
function distances(M::Matrix{Set{AminoAcid}}, distinctAAsets::Set{Tuple{Set{AminoAcid}, Set{AminoAcid}}})
	@info "$(loggingtime())\tdistances - main function" size(M) = size(M) distinctAAsets
	nrows = size(M, 1)
	Δmax = size(M, 2)  # maximal possible difference between any two proteins
	uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
	# for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
	Δs = tcollect(i_js_distances(M, distinctAAsets, uint, i) for i ∈ 1:nrows)
	# merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
	# Δ = SymmetricPacked(reduce(vcat, Δs))
	Δ = SymmetricPacked(reduce(vcat, Δs), :L)  # define the matrix by supplying its lower triangular part
	return Δ
end


function distances(M::Matrix{Set{AminoAcid}}, ::Empty{Set})
	@info "$(loggingtime())\tdistances - main function - empty distinctAAsets" size(M) = size(M)
	# the number of proteins to compare agaisnt each other
	nrows = size(M, 1)
	# since the are no distinctAAsets (distinctAAsets::Empty{Set}),
	# the distance between each two proteins is zero
	Δ = SymmetricPacked(zeros(UInts[1], nrows, nrows))
	return Δ
end


# example_unchosen_prot_df[!, findfirst(x->x=="254:257(K)", names(example_unchosen_prot_df)):end]
# bad_chosen_df[!, findfirst(x->x=="254:257(K)", names(bad_chosen_df)):end]

# # first row is our example unchosen protein, 
# # and the rest are all chosen proteins which we expect that
# # most of them will have distance == 0 with it
# M = Matrix(vcat(
# 	example_unchosen_prot_df[!, findfirst(x->x=="254:257(K)", names(example_unchosen_prot_df)):end],
# 	bad_chosen_df[!, findfirst(x->x=="254:257(K)", names(bad_chosen_df)):end]
# ))

# r1 = M[1, :]
# r2 = M[2, :]

# example_Δ = distances(M)


# all(.~isdisjoint.(r1, r2))


# disjoint_rows_indices = [
# 	i
# 	for (i, r2) in enumerate(eachrow(M[2:end, :]))
# 	if all(.~isdisjoint.(r1, r2))
# ]




# @assert isdisjoint.(
# 	[[1, 2], [2, 3, 4]], 
# 	[[2, 3, 4], [5, 6]]
# 	)  == [0, 1]
# @assert isdisjoint.(
# 	[[1, 2], [2, 3, 4]], 
# 	[[5], [5, 6]]
# 	) == [1, 1]





# for (s1, s2) in distinctAAsets
# 	@assert (s2, s1) in distinctAAsets
# end

# d1 = i_js_distances(M, distinctAAsets, uint, 1)

# all(length.(Δs) .== collect(range(nrows,1;step=-1)))

# temp_Δs = [
# 	[0, 1, 2],
# 	[0, 4],
# 	[0]
# ]

# reduce(vcat, temp_Δs)

# SymmetricPacked(reduce(vcat, temp_Δs))

# SymmetricPacked(reduce(vcat, temp_Δs), :L)



"""
	i_js_distances(M, distinctAAsets, uint, i)

Return `δᵢ` which is a vector of distances between `rowᵢ` to every `rowⱼ` in `M`, where `i <= j <= size(M, 1)`.  
A distance between `rowᵢ` to a `rowⱼ` is determined by to the number of corresponding positions in which the 
two have distinct sets of possible amino acid (naturally, the distance between `rowᵢ` to itself is 0).  
All possible distinct sets  between any two cells are given at `distinctAAsets`.  
`uint` is the type of unsigned int of the distances.
"""
function i_js_distances(
	M::Matrix{Set{AminoAcid}}, distinctAAsets, uint, i,
)
	@debug "$(loggingtime())\ti_js_distances - start" size(M) = size(M) i

	@views rowᵢ = M[i, :]
	nrows = size(M, 1)
	δᵢ = Vector{uint}(undef, nrows - i + 1)
	# the distance between `rowᵢ` to itself
	δᵢ[1] = 0
	# i == nrows && return δᵢ
	if i == nrows
		@debug "$(loggingtime())\ti_js_distances - end" size(M) = size(M) i
		return δᵢ
	end
	# the distances between `rowᵢ` to every `rowⱼ` where `i < j`
	for j ∈ i+1:nrows   # todo use @inbounds? https://docs.julialang.org/en/v1/base/base/#Base.@inbounds
		@views rowⱼ = M[j, :]
		δᵢⱼ::uint = 0   # distance between `rowᵢ` and `rowⱼ`
		for (x, y) ∈ zip(rowᵢ, rowⱼ)
			# increase `δᵢⱼ` for each position in which
			# the two rows of proteins have distinct sets of possible amino acids
			if (x, y) ∈ distinctAAsets
				δᵢⱼ += 1
			end
		end
		δᵢ[j-i+1] = δᵢⱼ
	end
	@debug "$(loggingtime())\ti_js_distances - end" size(M) = size(M) i
	return δᵢ
end


"""
	allargmins(A)

Return a vector filled with each index `i` such that `A[i] == minimum(A)`.
"""
function allargmins(A)

	# todo As of now, A must be a vector and so `1 <= i <= length(A)`
	# todo wrap into a standalone package

	minval = A[1]

	uint = smallestuint(length(A))

	argmins = Vector{uint}(undef, 1)
	argmins[1] = 1

	for i ∈ eachindex(A[2:end])
		j = i + 1
		if A[j] == minval
			# additional argmin
			push!(argmins, j)
		elseif A[j] < minval
			# new argmin
			minval = A[j]
			argmins = Vector{uint}(undef, 1)
			argmins[1] = j
		end
	end

	return argmins
end


"Get indices of eligible solutions from `distinctdf` on which expression levels should be calculated."
function choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)
	@info "$(loggingtime())\tchoosesolutions" fractions algs onlymaxdistinct

	_distinctdf = distinctdf[in.(distinctdf.Algorithm, Ref(algs)).&in.(distinctdf.Fraction, Ref(fractions)), :]

	if onlymaxdistinct
		maxdistinct = maximum(_distinctdf[!, "NumUniqueSamples"])
		_distinctdf = subset(
			distinctdf,
			"NumUniqueSamples" => x -> x .== maxdistinct,
		)

	end

	solutions = _distinctdf[:, "Index"]
end


colhasambiguousAAs(col) = any(length.(col) .> 1)


"""
	Calculate the entropy of the disfferent AAs distribution in a position.
	`col` is a vector of sets of AminoAcids, where each set represents the possible AAs in the position.
"""
function colAAentropy(col::Vector{Set{AminoAcid}})
	# keep only unambiguous AAs for this col
	col = col[length.(col).==1]
	# counts = collect(values(counter(col)))
	counts = values(counter(col))
	# S = - sum(p_i * log2(p_i))
	# where p_i = count_i / sum(counts)
	totalcount = sum(counts)
	if totalcount == 0
		return 0.0
	end
	p = counts ./ totalcount
	return -sum(p .* log2.(p))
end


"""
	Calculate the entropy of the amino acids in a protein.
	`protAAs` is a vector of sets of AminoAcids, where each set represents the possible AAs in the position.
	The entropy is the sum of the entropies of each unambiguous position of the protein.
"""
protentropy(protAAs::Vector{Set{AminoAcid}}, s_cols::Vector{Float64}) = sum(
	[
	s_col
	for (s_col, AAs) in zip(s_cols, protAAs)
	if length(AAs) == 1
]
)


"""
	Annotate proteins with sufficient entropy.
	`allprotsdf` is a DataFrame of proteins, where the first column is the protein name and the last ones are the amino acids in each position.
	`firstcolpos` is the position of the first column with amino acids.
	Returns a DataFrame with an additional column "SufficientEntropy" indicating whether the protein has 
	sufficient entropy, and the updated `firstcolpos` position.

	How do we determine whether a protein has sufficient entropy?
	Assuming we have N ambiguous positions (i.e., positions where for at least one protein, we have two possible AAs),
	each having 2/3 unambiguous options, with known probabilities.
	For each of these positions, the entropy is -sum(p * log2(p)), 
	where p is the probability of each unambiguous AA in the position.
	The total entropy (approximation, neglecting correlations) s_0 is the sum of all S values for all positions.
	For each protein, its entropy s_prot is the sum of the S values over the AA for which it has known values. 
	If this sum s_prot > S_0 / 2, this protein may be reassigned.
"""
function findprotswithsufficiententropy!(
	allprotsdf, firstcolpos, samplename,
)
	@info "$(loggingtime())\tfindprotswithsufficiententropy" samplename firstcolpos

	df = allprotsdf[:, firstcolpos:end]
	df = df[
		:,
		[i for (i, col) in enumerate(eachcol(df)) if colhasambiguousAAs(col)],
	]

	s_cols = colAAentropy.(eachcol(df))
	s_0 = sum(s_cols)
	s_prots = protentropy.(Vector.(eachrow(df)), Ref(s_cols))
	sufficiententropy = s_prots .> s_0 / 2  # the proteins that have enough entropy to be accepted

	proteincolpos = findfirst(names(allprotsdf) .== "Protein")
	insertcols!(allprotsdf, proteincolpos + 1, "SufficientEntropy" => sufficiententropy)
	firstcolpos += 1

	return allprotsdf, firstcolpos
end


possibly_empty_array_sum(arr) = isempty(arr) ? 0 : sum(arr)


function one_solution_additional_assignment_considering_available_reads(
	distinctdf, allprotsdf, firstcolpos, Δ, solution, readsdf, samplename,
	considerentropy::Bool = false,
)
	@info "$(loggingtime())\tone_solution_additional_assignment_considering_available_reads" samplename solution considerentropy

	solutionrow = distinctdf[solution, :]

	baseallprotsdf = deepcopy(allprotsdf[:, begin:firstcolpos-1])

	# if "AvailableReads" ∈ names(solutionrow) && solutionrow["Fraction"] < 1.0
	if "AvailableReads" ∈ names(solutionrow)
		availablereads = solutionrow["AvailableReads"]
		allreadsperprotein = baseallprotsdf[!, "Reads"]
		availablereadsperprotein = [
			[read for read ∈ reads if read ∈ availablereads]
			for reads ∈ allreadsperprotein
		]
		baseallprotsdf[:, "Reads"] .= availablereadsperprotein
		baseallprotsdf[:, "NumOfReads"] .= length.(availablereadsperprotein)
		# filter out proteins not supported by currently-available reads
		baseallprotsdf = filter("NumOfReads" => x -> x > 0, baseallprotsdf)
	end

	# calculate the mean na/edited/unedited positions per each protein's reads - considering the available reads

	cols_to_avg_in_reads_df = ["EditedPositions", "UneditedPositions", "NAPositions"]

	cols_to_avg_in_prots_df = ["EditedPositions", "UneditedPositions"]
	"NAPositions" in names(baseallprotsdf) ? push!(cols_to_avg_in_prots_df, "NAPositions") : push!(cols_to_avg_in_prots_df, "MeanNAPositions")

	for (col_to_avg_in_reads_df, col_to_avg_in_prots_df) in zip(cols_to_avg_in_reads_df, cols_to_avg_in_prots_df)
		mean_col_measured_per_proteins_reads = Vector{Float64}(undef, nrow(baseallprotsdf))
		@threads for idx in 1:nrow(baseallprotsdf)
			protein_row = baseallprotsdf[idx, :]
			protein_reads = protein_row["Reads"]
			col_measured_per_protein_reads = readsdf[readsdf[!, "Read"].∈Ref(protein_reads), col_to_avg_in_reads_df]
			mean_col_measured_per_protein_reads = mean(col_measured_per_protein_reads)
			mean_col_measured_per_proteins_reads[idx] = mean_col_measured_per_protein_reads
		end
		baseallprotsdf[!, col_to_avg_in_prots_df] .= mean_col_measured_per_proteins_reads
		mean_col_to_avg_in_prots_df = occursin("Mean", col_to_avg_in_prots_df) ? col_to_avg_in_prots_df : "Mean$col_to_avg_in_prots_df"
		mean_col_to_avg_in_prots_df in names(baseallprotsdf) || rename!(baseallprotsdf, col_to_avg_in_prots_df => mean_col_to_avg_in_prots_df)
	end

	prots_in_solution = solutionrow["UniqueSamples"]

	chosendf = filter("Protein" => x -> x ∈ prots_in_solution, baseallprotsdf)
	unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, baseallprotsdf)

	# counter(chosendf[:, "SufficientEntropy"])
	# counter(unchosendf[:, "SufficientEntropy"])

	if considerentropy
		# filter out unchosen proteins with insufficient entropy
		unchosendf = filter("SufficientEntropy" => x -> x, unchosendf)
	end

	chosenindices = chosendf[:, "Index"] # indices of chosen proteins in the complete Δ matrix
	unchosenindices = unchosendf[:, "Index"] # indices of unchosen proteins in the complete Δ matrix

	# before continuing any further,
	# validate that the distance between any two chosen proeins (two distinct proteins)
	# is at least 1 
	# (assuming there are at least 2 chosen ones)
	if length(chosenindices) >= 2
		chosenprot_chosenprot_distances = Δ[chosenindices, chosenindices] # the distances between the chosen proteins to themselves
		chosenprot_minimum_distances = minimum.(
			vcat(row[begin:i-1], row[i+1:end])
			for (i, row) in enumerate(eachrow(chosenprot_chosenprot_distances))
		) # the smallest distance between each chosen prot to all others
		@assert all(chosenprot_minimum_distances .>= 1)
	end

	# also report unchosen prots which could have been included in the MIS as well
	unchosenprot_chosenprot_distances = Δ[unchosenindices, chosenindices] # the distances between the unchosen proteins to chosen proteins
	unchosenprot_chosenprot_minimum_distances = minimum.(eachrow(unchosenprot_chosenprot_distances))
	# these are the candidates - unchosen prots with distance > 0 to each chosen prot
	unchosen_prots_with_min_d_1_rel_indices = findall(unchosenprot_chosenprot_minimum_distances .!== 0x00)

	if isempty(unchosen_prots_with_min_d_1_rel_indices)
		missing_chosen_prots = []
	else
		unchosen_prots_with_min_d_1_abs_indices = unchosenindices[unchosen_prots_with_min_d_1_rel_indices]

		# baseallprotsdf[unchosen_prots_with_min_d_1_abs_indices, "Protein"]


		# chosen_prots_with_min_d_to_missing_unchosen_rel_indices = []
		# for row in eachrow(Δ[unchosen_prots_with_min_d_1_abs_indices, chosenindices])
		# 	row_min = minimum(row)
		# 	push!(
		# 		chosen_prots_with_min_d_to_missing_unchosen_rel_indices,
		# 		findall(row .== row_min)
		# 	)
		# end
		# describe(length.(chosen_prots_with_min_d_to_missing_unchosen_rel_indices))

		# f = Figure(size = (400, 400))
		# xs = length.(chosen_prots_with_min_d_to_missing_unchosen_rel_indices)
		# ax = Axis(f[1, 1],
		# 	xlabel = "Closeset chosen proteins",
		# 	ylabel = "Missing proteins",
		# 	# yscale = log10,
		# 	# title = "Reassignments per unchosen protein",
		# 	limits = (0, nothing, 0, nothing),
		# 	xticks = 0:20:maximum(xs)+20,
		# 	# yticks = 0:50:600,
		# )
		# hist!(ax, xs, color = :red, strokecolor = :black, strokewidth = 1)
		# f

		# chosen_prots_with_min_d_to_missing_unchosen_rel_indices_counter = counter(
		# 	vcat(chosen_prots_with_min_d_to_missing_unchosen_rel_indices...)
		# )
		# df = DataFrame(ChosenProtRelIndex = vcat(chosen_prots_with_min_d_to_missing_unchosen_rel_indices...))
		# df = combine(
		# 	groupby(df, "ChosenProtRelIndex"),
		# 	nrow
		# )

		# f = Figure(size = (400, 400))
		# xs = df[!, :nrow]
		# ax = Axis(f[1, 1],
		# 	xlabel = "Missing proteins being closest\nto a chosen protein",
		# 	ylabel = "Chosen proteins",
		# 	# yscale = log10,
		# 	# title = "Reassignments per unchosen protein",
		# 	# limits = (0, nothing, 0, nothing),
		# 	xticks = 1:maximum(xs),
		# 	yticks = 0:50:600,
		# )
		# hist!(ax, xs, color = :red, strokecolor = :black, strokewidth = 1)
		# f

		final_unchosen_prots_with_min_d_1_abs_indices = unchosen_prots_with_min_d_1_abs_indices[
			sum.(eachrow(Δ[unchosen_prots_with_min_d_1_abs_indices, unchosen_prots_with_min_d_1_abs_indices])).==0
		]
		# if final_unchosen_prots_with_min_d_1_abs_indices is empty,
		# it means that all unchosen proteins are distinct from oneanother,
		# so only one of them can be added to the chsoen proteins' set - pick it at random
		if isempty(final_unchosen_prots_with_min_d_1_abs_indices)
			push!(
				final_unchosen_prots_with_min_d_1_abs_indices,
				only(sample(unchosen_prots_with_min_d_1_abs_indices, 1)),
			)
		end
		# missing_chosen_prots = baseallprotsdf[final_unchosen_prots_with_min_d_1_abs_indices, "Protein"]
		missing_chosen_prots = subset(baseallprotsdf, "Index" => ByRow(x -> x ∈ final_unchosen_prots_with_min_d_1_abs_indices))[!, "Protein"]

	end
	newsolutionrow = DataFrame(solutionrow)
	newsolutionrow[:, "MissingUniqueSamples"] = [missing_chosen_prots]
	newsolutionrow[:, "NumMissingUniqueSamples"] = length.(newsolutionrow[:, "MissingUniqueSamples"])

	insertcols!(
		chosendf,
		3,
		"#Solution" => solutionrow["Index"],
		"Fraction" => solutionrow["Fraction"],
		"FractionRepetition" => solutionrow["FractionRepetition"],
		"Algorithm" => solutionrow["Algorithm"],
		"AlgorithmRepetition" => solutionrow["AlgorithmRepetition"],
	)

	for unchosenprot ∈ eachrow(unchosendf)  # a row of unchosen protein relative to the chosen proteins in the solution

		# unchosenprot = unchosendf[1, :]  # TODO - comment out - an example row of unchosen protein
		# unchosenprot = unchosendf[2, :]  

		unchosenprot_index = unchosenprot["Index"] # the row number of the unchosen protein in the distances matrix Δ

		unchosenprot_distances = Δ[unchosenprot_index, chosenindices] # the distances between the unchosen protein `unchosenprot` to all chosen proteins

		unchosenprot_distances_argmins = allargmins(unchosenprot_distances) # indices of minimum distances

		minchosendf = chosendf[unchosenprot_distances_argmins, :] # chosen proteins with minimum distance to the unchosen protein `unchosenprot`

		unchosenreads = unchosenprot["NumOfReads"] # the number of reads `unchosenprot` divides between the chosen proteins that are closest to it

		equal_addition = unchosenreads / size(minchosendf, 1)
		weighted_additions = unchosenreads .* minchosendf[:, "NumOfReads"] ./ sum(minchosendf[:, "NumOfReads"])

		@assert isapprox(equal_addition * size(minchosendf, 1), sum(weighted_additions))

		chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"] .+= equal_addition
		chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"] .+= weighted_additions

		newsupportingreads = unchosenprot["Reads"]
		newsupportingprotein = unchosenprot["Protein"]

		if !isassigned(newsupportingreads)
			println("solution: ", solution)
			println("newsupportingreads: ", newsupportingreads)
			println("newsupportingprotein: ", newsupportingprotein)
			# throw(ErrorException("newsupportingreads or newsupportingprotein is not assigned"))
			break
		end

		for existingsupportingreads ∈ chosendf[unchosenprot_distances_argmins, "AdditionalSupportingReadsIDs"]
			push!(existingsupportingreads, newsupportingreads)
		end

		@assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution, Protein: $(unchosenprot["Protein"])"""
		@assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution, Protein: $(unchosenprot["Protein"])"""

		for existingsupportingproteins ∈ chosendf[unchosenprot_distances_argmins, "AdditionalSupportingProteinsIDs"]
			push!(existingsupportingproteins, newsupportingprotein)
		end

		@assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution, Protein: $(unchosenprot["Protein"])"""
		@assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution, Protein: $(unchosenprot["Protein"])"""


		min_distance = minimum(unchosenprot_distances)
		# chosendf[unchosenprot_distances_argmins, "AdditionalSupportingProteinsDistances"] .= [[] for i in eachindex(unchosenprot_distances_argmins)]
		push!.(chosendf[unchosenprot_distances_argmins, "AdditionalSupportingProteinsDistances"], min_distance)

		mean_na_positions = unchosenprot["MeanNAPositions"]
		push!.(chosendf[unchosenprot_distances_argmins, "AdditionalSupportingProteinsMeanNAPositions"], mean_na_positions)

		push!.(chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReadsContributionPerProtein"], equal_addition)

		for (existingsupportingreadscontrib, weighted_addition) ∈ zip(
			chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReadsContributionPerProtein"],
			weighted_additions,
		)
			push!(existingsupportingreadscontrib, weighted_addition)
		end

	end


	@assert all(
		isapprox(
			possibly_empty_array_sum.(chosendf[:, "AdditionalEqualSupportingReadsContributionPerProtein"]),
			chosendf[:, "AdditionalEqualSupportingReads"],
		),
	)
	@assert all(
		isapprox(
			possibly_empty_array_sum.(chosendf[:, "AdditionalWeightedSupportingReadsContributionPerProtein"]),
			chosendf[:, "AdditionalWeightedSupportingReads"],
		),
	)

	chosendf[:, "TotalEqualSupportingReads"] .= chosendf[:, "NumOfReads"] .+ chosendf[:, "AdditionalEqualSupportingReads"]
	chosendf[:, "TotalWeightedSupportingReads"] .= chosendf[:, "NumOfReads"] .+ chosendf[:, "AdditionalWeightedSupportingReads"]
	chosendf[:, "AdditionalSupportingProteins"] .= length.(chosendf[:, "AdditionalSupportingProteinsIDs"])

	sort!(chosendf, "AdditionalSupportingProteins")

	transform!(chosendf, :AdditionalSupportingProteinsIDs => ByRow(x -> join(x, ",")) => :AdditionalSupportingProteinsIDs)
	transform!(chosendf, :AdditionalSupportingReadsIDs => ByRow(x -> join(join.(x, ","), ";")) => :AdditionalSupportingReadsIDs)

	# chosendf[:, :AdditionalSupportingProteinsIDs]
	# transform(chosendf, :AdditionalSupportingReadsIDs => ByRow(x -> join(join.(x, ","), ";")) => :AdditionalSupportingReadsIDs)[:, :AdditionalSupportingReadsIDs]

	@assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution"""
	@assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution"""

	# return chosendf
	return chosendf, newsolutionrow

end


function additional_assignments(
	distinctdf, allprotsdf, firstcolpos, Δ, solutions, readsdf, samplename, considerentropy,
)
	@info "$(loggingtime())\tadditional_assignments" samplename solutions considerentropy
	results = map(solutions) do solution
		try
			one_solution_additional_assignment_considering_available_reads(
				distinctdf, allprotsdf, firstcolpos, Δ, solution, readsdf, samplename,
				considerentropy,
			)
		catch e
			@warn "Error in additional_assignments:" e solution
			missing
		end
	end
	return results
end


function run_sample(
	distinctfile, allprotsfile, samplename, postfix_to_add,
	firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
	maxmainthreads, outdir, algs, onlymaxdistinct,
	readsfile,
	substitutionmatrix::Union{SubstitutionMatrix, Nothing},
	similarityscorecutoff::Int64,
	similarityvalidator::Function,
	aagroups::Union{Dict{AminoAcid, String}, Nothing},
	considerentropy::Bool,
)
	@info "$(loggingtime())\trun_sample" distinctfile allprotsfile readsfile samplename considerentropy

	distinctdf = prepare_distinctdf(
		distinctfile, delim, innerdelim, truestrings, falsestrings,
	)

	readsdf = prepare_readsdf(readsfile, delim)

	allprotsdf, firstcolpos = prepare_allprotsdf!(
		allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos,
		readsdf,
	)

	if considerentropy
		allprotsdf, firstcolpos = findprotswithsufficiententropy!(
			allprotsdf, firstcolpos, samplename,
		)
	end

	# the possible amino acids each protein has in each position
	M = Matrix(allprotsdf[:, firstcolpos:end])
	# the distances between any two proteins according to `M`
	Δ = begin
		if substitutionmatrix !== nothing
			distances(M, substitutionmatrix, similarityscorecutoff, similarityvalidator)
		elseif aagroups !== nothing
			distances(M, aagroups)
		else
			distances(M)
		end
	end

	# considering only desired solutions (rows' indices)
	solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)

	# allsubsolutions = collect(Iterators.partition(solutions, maxmainthreads))

	# minmainthreads = minimum([Int(Threads.nthreads() / 5), Int(length(solutions) / 4)])
	minmainthreads = max(
		1,
		min(
			Int(round(Threads.nthreads() / 5)),
			Int(round(length(solutions) / 4)),
		),
	)
	allsubsolutions = collect(Iterators.partition(solutions, minmainthreads))

	results = tcollect(
		additional_assignments(
			distinctdf, allprotsdf, firstcolpos, Δ, subsolutions, readsdf, samplename,
			considerentropy,
		)
		for subsolutions ∈ allsubsolutions
	)
	# finalresults = vcat(Iterators.flatten(results)...)
	# finalresults = vcat(skipmissing(Iterators.flatten(results)...))
	finalresults = vcat((skipmissing(Iterators.flatten(results))...))

	# save these into seperate files, 
	# the first being `outfile` and the second an "updated" version of the `distinctfile`
	finalexpresults = vcat([result[1] for result in finalresults]...)
	newdistinctdf = vcat([result[2] for result in finalresults]...)


	# save the expression results
	finalexpresults[!, "Reads"] .= join.(finalexpresults[!, "Reads"], innerdelim)
	roundelements(x, digits = 5) = round.(x; digits)
	stringifyelements(x) = string.(x)
	intifyelements(x) = Int.(x)
	transform!(finalexpresults, :AdditionalEqualSupportingReadsContributionPerProtein => x -> stringifyelements.(roundelements.(x)), renamecols = false)
	finalexpresults[!, "AdditionalEqualSupportingReadsContributionPerProtein"] .= join.(finalexpresults[!, "AdditionalEqualSupportingReadsContributionPerProtein"], innerdelim)
	transform!(finalexpresults, :AdditionalWeightedSupportingReadsContributionPerProtein => x -> stringifyelements.(roundelements.(x)), renamecols = false)
	finalexpresults[!, "AdditionalWeightedSupportingReadsContributionPerProtein"] .= join.(finalexpresults[!, "AdditionalWeightedSupportingReadsContributionPerProtein"], innerdelim)
	transform!(finalexpresults, :AdditionalSupportingProteinsDistances => x -> stringifyelements.(intifyelements.(x)), renamecols = false)
	finalexpresults[!, "AdditionalSupportingProteinsDistances"] .= join.(finalexpresults[!, "AdditionalSupportingProteinsDistances"], innerdelim)
	transform!(finalexpresults, :AdditionalSupportingProteinsMeanNAPositions => x -> stringifyelements.(roundelements.(x)), renamecols = false)
	finalexpresults[!, "AdditionalSupportingProteinsMeanNAPositions"] .= join.(finalexpresults[!, "AdditionalSupportingProteinsMeanNAPositions"], innerdelim)

	expoutfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels$postfix_to_add.csv")
	CSV.write(expoutfile, finalexpresults; delim)

	# save the updated distinctdf
	empty_arr_to_empty_string(arr) = isempty(arr) ? "" : arr
	newdistinctdf[!, "MissingUniqueSamples"] .= join.(empty_arr_to_empty_string.(newdistinctdf[!, "MissingUniqueSamples"]))
	newdistinctdf = outerjoin(distinctdf, newdistinctdf, on = names(distinctdf))
	sort!(newdistinctdf, "Index")
	# newdistinctdf[!, "MissingUniqueSamples"] .= ifelse.(ismissing.(newdistinctdf[!, "MissingUniqueSamples"]), [[]], newdistinctdf[!, "MissingUniqueSamples"])
	# newdistinctdf[!, "NumMissingUniqueSamples"] .= ifelse.(ismissing.(newdistinctdf[!, "NumMissingUniqueSamples"]), 0, newdistinctdf[!, "NumMissingUniqueSamples"])
	newdistinctdf[!, "UniqueSamples"] .= join.(newdistinctdf[!, "UniqueSamples"], innerdelim)
	newdistinctdf[!, "AvailableReads"] .= join.(newdistinctdf[!, "AvailableReads"], innerdelim)
	distinctfilepostfixstart = findlast('.', distinctfile)
	updateddistinctfile = distinctfile[begin:distinctfilepostfixstart-1] * ".Updated" * "$postfix_to_add" * distinctfile[distinctfilepostfixstart:end]
	CSV.write(updateddistinctfile, newdistinctdf; delim)
end


function main(
	distinctfiles, allprotsfiles, samplenames,
	postfix_to_add,
	firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
	maxmainthreads, outdir, algs, onlymaxdistinct,
	gcp, shutdowngcp,
	readsfiles,
	substitutionmatrix::Union{SubstitutionMatrix, Nothing},
	similarityscorecutoff::Int64,
	similarityvalidator::Function,
	aagroups::Union{Dict{AminoAcid, String}, Nothing},
	considerentropy::Bool,
	logtostdout::Bool,
	minloglevel,
)
	# logfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expression.logger.log"

	# logfile !== nothing && global_logger(FileLogger(logfile; append = true))

	# logfile !== nothing && global_logger(MinLevelLogger(FileLogger(logfile; append = true), Logging.Info))
	# if logfile === nothing
	# 	logfile = "$outdir/expressionlevels$postfix_to_add.$(Dates.format(now(), "dd.mm.YYYY")).log"
	# end
	if logtostdout
		global_logger(MinLevelLogger(ConsoleLogger(stdout, minloglevel), minloglevel))
	else
		logfile = "$outdir/expressionlevels$postfix_to_add.$(Dates.format(now(), "dd.mm.YYYY")).log"
		# global_logger(MinLevelLogger(FileLogger(logfile; append = true), minloglevel))
		global_logger(MinLevelLogger(FileLogger(logfile), minloglevel))
	end


	@info "$(loggingtime())\tmain" distinctfiles allprotsfiles readsfiles samplenames postfix_to_add firstcolpos delim innerdelim truestrings falsestrings fractions maxmainthreads outdir algs onlymaxdistinct gcp shutdowngcp substitutionmatrix similarityscorecutoff similarityvalidator aagroups considerentropy logtostdout minloglevel

	length(distinctfiles) == length(allprotsfiles) == length(readsfiles) == length(samplenames) || error("Unequal input files' lengths!")

	# run each sample using the `run_sample` function
	for (distinctfile, allprotsfile, samplename, readsfile) ∈ zip(distinctfiles, allprotsfiles, samplenames, readsfiles)
		run_sample(
			distinctfile, allprotsfile, samplename, postfix_to_add,
			firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
			maxmainthreads, outdir, algs, onlymaxdistinct,
			readsfile,
			substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups,
			considerentropy,
		)
	end

	@info "$(loggingtime())\tmain - finished sucessfully!"

	# shutdown gcp vm (if needed)
	gcp && shutdowngcp && run(`sudo shutdown`) # https://cloud.google.com/compute/docs/shutdownscript
end


"""
Define command-line arguments.
"""
function parsecmd()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--distinctfiles"
		help = "One or more csv files representing distinct unique proteins."
		nargs = '+'
		action = :store_arg
		# required = true
		# default = []
		"--allprotsfiles"
		help = "corresponding csv files representing unique proteins (w.r.t. `distinctfiles`)."
		nargs = '+'
		action = :store_arg
		"--allreadsfiles"
		help = "corresponding csv files representing the basic reads (w.r.t. `distinctfiles`)."
		nargs = '+'
		action = :store_arg
		# required = true
		# default = []
		"--samplenames"
		nargs = '+'
		action = :store_arg
		# required = true
		# default = []

		"--distinctfilesfofn"
		help = "One or more csv files representing distinct unique proteins."
		# nargs = '+'
		# action = :store_arg
		# required = true
		default = ""
		"--allprotsfilesfofn"
		help = "corresponding csv files representing unique proteins (w.r.t. `distinctfiles`)."
		# nargs = '+'
		# action = :store_arg
		# required = true
		default = ""
		"--allreadsfilesfofn"
		help = "corresponding csv files representing reads (w.r.t. `distinctfiles`)."
		# nargs = '+'
		# action = :store_arg
		# required = true
		default = ""
		"--samplenamesfile"
		# nargs = '+'
		# action = :store_arg
		# required = true
		default = ""

		"--postfix_to_add"
		help = "Add `postfix` to output files' names, e.g., `\$sample.DistinctUnique{Reads,Proteins}\$postfix.\$time.csv`."
		default = ""
		"--firstcolpos"
		help = "Int location of the first editing position column of each file in `allprotsfiles`."
		arg_type = Int
		default = 15
		"--delim"
		help = "Delimiter for input/output csv files."
		default = "\t"
		"--innerdelim"
		help = "Inner delimiter for cells with multiple values in input/output csv files."
		default = ","

		"--truestrings"
		help = "Possible representations of truestrings values in csv files. Used for `CSV.FILE`."
		nargs = '+'
		action = :store_arg
		default = ["TRUE", "True", "true"]
		"--falsestrings"
		help = "Possible representations of false values in csv files. Same for `CSV.FILE`."
		nargs = '+'
		action = :store_arg
		default = ["FALSE", "False", "false"]

		"--fractions"
		help = "Consider only solutions of desired fractions."
		default = [1.0]
		nargs = '+'
		action = :store_arg
		arg_type = Float64
		"--algs"
		help = "Consider only solutions achieved by desired algorithms."
		default = ["Ascending", "Descending"]
		nargs = '+'
		range_tester = x -> x ∈ ["Ascending", "Descending", "Unordered"]
		action = :store_arg
		"--onlymaxdistinct"
		help = "Consider only solution with maximum number of distinct proteins, within allowed fractions and algorithms."
		action = :store_true

		"--maxmainthreads"
		default = 30
		arg_type = Int

		"--outdir"
		help = "Write output files to this directory."
		required = true

		"--gcp"
		help = "Program is run on a google cloud VM."
		action = :store_true
		"--shutdowngcp"
		help = "Shutdown google cloud VM when the program ends."
		action = :store_true

		"--substitutionmatrix"
		help = "Use this substitution matrix as a stricter criteria for determination of distinct AAs. Use in conjuction with `datatype == Proteins`. Not compatible with `aagroups`. Use any matrix in https://github.com/KaparaNewbie/BioAlignments.jl/blob/master/src/submat.jl."
		arg_type = Symbol
		"--similarityscorecutoff"
		help = "See `similarityvalidator` below."
		arg_type = Int
		default = 0
		"--similarityvalidator"
		help = "Use this opeartor to determine similarty of AA change, e.g., whether `5 >= similarityscorecutoff`."
		arg_type = Symbol
		default = :(>=)
		"--aagroups"
		help = "Use predifined AAs classification as a stricter criteria for determination of distinct AAs. Use in conjuction with `datatype == Proteins`. Not compatible with `substitutionmatrix`."
		arg_type = Symbol
		range_tester = x -> x ∈ [:AA_groups, :AA_groups_Miyata1979]
		"--considerentropy"
		help = "If `true`, will only reassign unchosen proteins with sufficient entropy."
		action = :store_true

		"--logtostdout"
		help = """
		Write log messages to `stdout` instead of a designated file 
		(its path will be configured by `outdir` and `postfix_to_add`)
		However, the time flushing to `stdout` is not guaranteed, so it is recommended 
		to let the program flush the logging to a log file to follow in real time.
		"""
		action = :store_true
		"--minloglevel"
		# default = Logging.Info
		# arg_type = LogLevel
		default = 0
		arg_type = Int
		help = "Minimum log level to be logged. Default is 0 (Info). Use -1000 for debugging."

	end
	return parse_args(s)
end


function CLI_main()
	# read command-line args
	parsedargs = parsecmd()

	distinctfiles = parsedargs["distinctfiles"]
	allprotsfiles = parsedargs["allprotsfiles"]
	allreadsfiles = parsedargs["allreadsfiles"]
	samplenames = parsedargs["samplenames"]

	distinctfilesfofn = parsedargs["distinctfilesfofn"]
	allprotsfilesfofn = parsedargs["allprotsfilesfofn"]
	allreadsfilesfofn = parsedargs["allreadsfilesfofn"]
	samplenamesfile = parsedargs["samplenamesfile"]

	postfix_to_add = parsedargs["postfix_to_add"]

	firstcolpos = parsedargs["firstcolpos"]
	delim = parsedargs["delim"]
	innerdelim = parsedargs["innerdelim"]
	truestrings = parsedargs["truestrings"]
	falsestrings = parsedargs["falsestrings"]

	fractions = parsedargs["fractions"]
	# println("fractions = $fractions, typeof(fractions) = $typeof(fractions)")
	# fractions = fractions isa Vector{Float64} ? fractions : parse.(Float64, fractions)
	if !(fractions isa Vector{Float64})
		fractions = Float64.(fractions)
	end
	algs = parsedargs["algs"]
	onlymaxdistinct = parsedargs["onlymaxdistinct"]

	maxmainthreads = parsedargs["maxmainthreads"]
	outdir = parsedargs["outdir"]

	gcp = parsedargs["gcp"]
	shutdowngcp = parsedargs["shutdowngcp"]

	substitutionmatrix = eval(parsedargs["substitutionmatrix"])
	similarityscorecutoff = parsedargs["similarityscorecutoff"]
	similarityvalidator = eval(parsedargs["similarityvalidator"])
	aagroups = eval(parsedargs["aagroups"])

	considerentropy = parsedargs["considerentropy"]

	logtostdout = parsedargs["logtostdout"]
	minloglevel = LogLevel(parsedargs["minloglevel"])

	if distinctfilesfofn != ""
		distinctfiles = split(readline(distinctfilesfofn), " ")
	end
	if allprotsfilesfofn != ""
		allprotsfiles = split(readline(allprotsfilesfofn), " ")
	end
	if allreadsfilesfofn != ""
		allreadsfiles = split(readline(allreadsfilesfofn), " ")
	end

	if samplenamesfile != ""
		samplenames = split(readline(samplenamesfile), " ")
	end

	if truestrings !== nothing
		truestrings = convert(Vector{String}, truestrings)
	end
	if falsestrings !== nothing
		falsestrings = convert(Vector{String}, falsestrings)
	end

	# run
	main(
		distinctfiles, allprotsfiles, samplenames,
		postfix_to_add,
		firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
		maxmainthreads, outdir, algs, onlymaxdistinct,
		gcp, shutdowngcp,
		allreadsfiles,
		substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups,
		considerentropy,
		logtostdout, minloglevel,
	)
end


if abspath(PROGRAM_FILE) == @__FILE__
	CLI_main()
end


# distinctfile =  "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/DistinctProteins/comp183648_c0_seq1.DistinctUniqueProteins.21.03.2024-22:37:27.csv"
# allprotsfile = "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/ProteinsFiles/comp183648_c0_seq1.unique_proteins.csv.gz"
# readsfile = "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq.Polished.Unclustered.TotalCoverage50.BQ30.AHL.BHAfterNoise.3/ReadsFiles/comp183648_c0_seq1.reads.csv.gz"
# samplename = "comp183648_c0_seq1"	
# firstcolpos = 16


# distinctfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.DistinctUniqueProteins.06.02.2024-09:29:20.csv"
# readsfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.reads.csv.gz"
# allprotsfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# samplename = "GRIA"
# firstcolpos = 15
# onlymaxdistinct = true
# considerentropy = true



# delim = "\t"
# innerdelim = ","
# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]
# onlymaxdistinct = true
# algs = ["Ascending", "Descending"]
# maxmainthreads = 30
# fractions = [1.0]

# substitutionmatrix = nothing
# aagroups = nothing
# similarityscorecutoff = 0
# similarityvalidator = :(>=)


# distinctdf = prepare_distinctdf(
# 	distinctfile, delim, innerdelim, truestrings, falsestrings,
# )

# readsdf = prepare_readsdf(readsfile, delim)

# allprotsdf, firstcolpos = prepare_allprotsdf!(
# 	allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos, readsdf,
# )

# if considerentropy
# 	allprotsdf, firstcolpos = findprotswithsufficiententropy!(
# 		allprotsdf, firstcolpos
# 	)
# end


# # the possible amino acids each protein has in each position
# M = Matrix(allprotsdf[:, firstcolpos:end])
# # the distances between any two proteins according to `M`
# Δ = begin
# 	if substitutionmatrix !== nothing
# 		distances(M, substitutionmatrix, similarityscorecutoff, similarityvalidator)
# 	elseif aagroups !== nothing
# 		distances(M, aagroups)
# 	else
# 		distances(M)
# 	end
# end



# # # # Int(maximum(Δ)) # maximum distance between any two proteins
# # # # Int(minimum(Δ)) # minimum distance between any two proteins


# solution = distinctdf[distinctdf[!, "NumUniqueSamples"].==maximum(distinctdf[!, "NumUniqueSamples"]), "Index"][1]
# @assert length(solution) == 1 # should be a single element - not a vector

# result = one_solution_additional_assignment_considering_available_reads(
# 				distinctdf, allprotsdf, firstcolpos, Δ, solution, readsdf, samplename,
# 				considerentropy
# 			)

# result = result[1]


# prots_in_solution = distinctdf[solution, "UniqueSamples"]

# chosendf = filter("Protein" => x -> x ∈ prots_in_solution, allprotsdf)
# unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, allprotsdf)

# unique_reassigned_prots = filter(x -> !isempty(x), unique(vcat(split.(result[!, "AdditionalSupportingProteinsIDs"], ",")...)))

# reassigned_prots_counter = counter(filter(x -> !isempty(x), vcat(split.(result[!, "AdditionalSupportingProteinsIDs"], ",")...)))


# reassigned_unchosendf = filter("Protein" => x -> x ∈ unique_reassigned_prots, unchosendf)

# unreassigned_unchosendf = filter("Protein" => x -> x ∉ unique_reassigned_prots, unchosendf)


# counter(reassigned_unchosendf[:, "SufficientEntropy"])
# counter(unreassigned_unchosendf[:, "SufficientEntropy"])


# # solution = 15

# considering only desired solutions (rows' indices)
# solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)
# solutions = choosesolutions(distinctdf, [0.2], algs, false)
# minmainthreads = max(
# 	1,
# 	min(
# 		Int(round(Threads.nthreads() / 5)),
# 		Int(round(length(solutions) / 4))
# 	)
# )
# # minmainthreads = 8
# allsubsolutions = collect(Iterators.partition(solutions, minmainthreads))





# logfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/expression.logger.log"
# logger = FileLogger(logfile; append = true)
# # logger = ConsoleLogger()
# global_logger(logger)
# results = tcollect(
# 	additional_assignments(distinctdf, allprotsdf, firstcolpos,
# 	Δ, subsolutions, readsdf, samplename)
# 	for subsolutions ∈ allsubsolutions
# )

# results = tcollect(
# 	additional_assignments(distinctdf, allprotsdf, firstcolpos,
# 	Δ, subsolutions, readsdf, samplename)
# 	for subsolutions ∈ allsubsolutions
# )


# # finalresults = vcat(Iterators.flatten(results)...)
# # finalresults = vcat(skipmissing(Iterators.flatten(results)...))
# finalresults = vcat((skipmissing(Iterators.flatten(results))...))

# finalresults[8][1]
# finalresults[8][2]

# size(finalresults)

# finalexpresults = vcat([result[1] for result in finalresults]...)

# newdistinctdf = vcat([result[2] for result in finalresults]...)
# sort(newdistinctdf[!,["Algorithm", "NumMissingUniqueSamples"]], "Algorithm")




# # plot distances between all protein pairs
# all_tri_ds_counter = counter(Δ.tri)
# diag_ds_counter = counter(diag(Δ))
# for (k, v) in diag_ds_counter
# 	all_tri_ds_counter[k] -= v
# end
# for (k, v) in all_tri_ds_counter
# 	if v == 0x00
# 		pop!(all_tri_ds_counter, k)
# 	end
# end
# distancesdf = DataFrame(Distance = Int.(keys(all_tri_ds_counter)), Count = collect(values(all_tri_ds_counter)))
# xs = distancesdf[!, "Distance"]
# ys = distancesdf[!, "Count"]
# f = Figure(size = (500, 500))
# ax = Axis(f[1, 1],
# 	# xlabel = "Distance", 
# 	xlabel = "Distance",
# 	ylabel = "Protein pairs",
# 	yscale = log10,
# 	title = "Distance between any two proteins",
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# barplot!(ax, xs, ys, color = :red, strokecolor = :black, strokewidth = 1)
# f





# # considering only desired solutions (rows' indices)
# solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)



# # allsubsolutions = collect(Iterators.partition(solutions, maxmainthreads))

# minmainthreads = minimum([Int(Threads.nthreads() / 5), Int(length(solutions) / 4)])
# allsubsolutions = collect(Iterators.partition(solutions, minmainthreads))


# # solution = 79
# solution = distinctdf[distinctdf[!, "NumUniqueSamples"].==maximum(distinctdf[!, "NumUniqueSamples"]), "Index"][1]
# @assert length(solution) == 1 # should be a single element - not a vector

# solutionrow = distinctdf[solution, :]

# # todo update these lines based on the updated function they are taken from

# baseallprotsdf = deepcopy(allprotsdf[:, begin:firstcolpos-1])

# # # if "AvailableReads" ∈ names(solutionrow) && solutionrow["Fraction"] < 1.0
# # if "AvailableReads" ∈ names(solutionrow)
# # 	availablereads = solutionrow["AvailableReads"]
# # 	allreadsperprotein = baseallprotsdf[!, "Reads"]
# # 	availablereadsperprotein = [
# # 		[read for read ∈ reads if read ∈ availablereads]
# # 		for reads ∈ allreadsperprotein
# # 	]
# # 	baseallprotsdf[:, "Reads"] .= availablereadsperprotein
# # 	baseallprotsdf[:, "NumOfReads"] .= length.(availablereadsperprotein)
# # 	baseallprotsdf = filter("NumOfReads" => x -> x > 0, baseallprotsdf)
# # end

# prots_in_solution = solutionrow["UniqueSamples"]

# chosendf = filter("Protein" => x -> x ∈ prots_in_solution, baseallprotsdf)
# unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, baseallprotsdf)

# chosenindices = chosendf[:, "Index"] # indices of chosen proteins in the complete Δ matrix
# unchosenindices = unchosendf[:, "Index"] # indices of unchosen proteins in the complete Δ matrix

# # prots_not_in_solution = setdiff(allprotsdf[!, "Protein"], prots_in_solution)




# # example_unchosen_prot = "27"
# # example_unchosen_prot_df = allprotsdf[allprotsdf[!, "Protein"] .== example_unchosen_prot, :]
# # example_unchosen_prot_index = unchosendf[unchosendf[!, "Protein"] .== example_unchosen_prot, "Index"]
# # example_unchosen_prot_distances = Δ[example_unchosen_prot_index, chosenindices] # the distances between the 

# # bad_chosen_indices = [i for (i, d) in zip(chosenindices, example_unchosen_prot_distances) if d !== 0x00]
# # bad_chosen_df = copy(allprotsdf[bad_chosen_indices, :])
# # insertcols!(
# # 	bad_chosen_df, 
# # 	3, 
# # 	"DistanceFrom$example_unchosen_prot" => example_unchosen_prot_distances[example_unchosen_prot_distances .!== 0x00]
# # 	)


# chosenprot_chosenprot_distances = Δ[chosenindices, chosenindices] # the distances between the chosen proteins to themselves
# chosenprot_minimum_distances = minimum.(
# 	vcat(row[begin:i-1], row[i+1:end])
# 	for (i, row) in enumerate(eachrow(chosenprot_chosenprot_distances))
# )
# @assert all(chosenprot_minimum_distances .>= 1)





# unchosenprot_chosenprot_distances = Δ[unchosenindices, chosenindices] # the distances between the unchosen proteins to chosen proteins

# minimum_distances = minimum.(eachrow(unchosenprot_chosenprot_distances))
# describe(minimum_distances)


# unchosen_prots_with_min_d_1_rel_indices = findall(minimum_distances .!== 0x00)
# unchosen_prots_with_min_d_1_abs_indices = unchosenindices[unchosen_prots_with_min_d_1_rel_indices]

# unchosen_prots_with_min_d_1_distances = Δ[unchosen_prots_with_min_d_1_abs_indices, chosenindices]
# describe(minimum.(eachrow(unchosen_prots_with_min_d_1_distances)))

# one_unchosen_prot_index = unchosen_prots_with_min_d_1_abs_indices[1]
# five_chosen_prots_indices = chosenindices[4:8]

# _ds = Int.(Δ[one_unchosen_prot_index, five_chosen_prots_indices])

# _df = allprotsdf[vcat(one_unchosen_prot_index, five_chosen_prots_indices), firstcolpos:end]



# # # plot minumum distance between each unchosen protein to the chosen proteins
# # f = Figure(size = (500, 500))
# # ax = Axis(f[1, 1],
# # 	# xlabel = "Distance", 
# # 	xlabel = "Min distance to a chosen protein",
# # 	ylabel = "Unchosen proteins",
# # 	yscale = log10,
# # 	title = "Min distance between unchosen proteins to chosen ones",
# # 	# limits = (0, nothing, 1, nothing)
# # 	# yticks = LogTicks(WilkinsonTicks(2))
# # 	# yticks = [10 ^ i for i in 1:9]
# # )
# # hist!(ax, minimum_distances, color = :red, strokecolor = :black, strokewidth = 1)
# # f

# unchosen_distances_counters = counter.(eachrow(unchosenprot_chosenprot_distances))

# per_unchosen_min_distance_counts = [
# 	unchosen_distances_counter[min_distance]
# 	for (unchosen_distances_counter, min_distance)
# 	in
# 	zip(unchosen_distances_counters, minimum_distances)
# ]
# describe(per_unchosen_min_distance_counts)

# # f = Figure(size = (500, 500))
# # ax = Axis(f[1, 1],
# # 	# xlabel = "Distance", 

# # 	xlabel = "Chosen proteins with min distance to unchosen protein",
# # 	ylabel = "Unchosen proteins",
# # 	yscale = log10,
# # 	title = "Reassignments per unchosen protein",
# # 	# limits = (0, nothing, 1, nothing)
# # 	# yticks = LogTicks(WilkinsonTicks(2))
# # 	# yticks = [10 ^ i for i in 1:9]
# # )
# # hist!(ax, per_unchosen_min_distance_counts, color = :red, strokecolor = :black, strokewidth = 1)
# # f


# # let's plot the following things:
# # 1 - how many different distances an unchosen protein has?
# # 2 - avg distance to chosen proteins
# # 3 - the 3 smallest distances between each unchosen protein to the chosen proteins
# # 4 - the number of chosen proteins per the 3 smallest distances for each unchosen protein



# # 1 - how many different distances an unchosen protein has?

# unique_distances_per_unchosen_prot = sort(length.(keys.(unchosen_distances_counters)))
# describe(unique_distances_per_unchosen_prot)

# # sum(unique_distances_per_unchosen_prot .== 2)
# # unchosen_distances_counters[unique_distances_per_unchosen_prot .== 2]

# unchosen_prots_unique_distances_counter = counter(unique_distances_per_unchosen_prot)
# sorted_unique_distances = sort(collect(keys(unchosen_prots_unique_distances_counter)))
# sorted_unique_distances_counts = [unchosen_prots_unique_distances_counter[d] for d in sorted_unique_distances]

# sorted_unique_distances_counts_cumsum = cumsum(sorted_unique_distances_counts)
# sorted_unique_distances_counts_cumsum_percent = 100 .* sorted_unique_distances_counts_cumsum ./ sorted_unique_distances_counts_cumsum[end]

# f = Figure(size = (500, 500))
# ax = Axis(f[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Num of unique distances to chosen proteins per unchosen protein",
# 	ylabel = "Unchosen proteins",
# 	yscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# hist!(ax, unique_distances_per_unchosen_prot, color = :red, strokecolor = :black, strokewidth = 1)
# f


# # f = Figure(size = (500, 500))
# # ax = Axis(f[1, 1],
# # 	# xlabel = "Distance", 

# # 	xlabel = "Num of unique distances to chosen proteins per unchosen protein",
# # 	ylabel = "Unchosen proteins (cumulative)",
# # 	# yscale = log10,
# # 	# title = "Reassignments per unchosen protein",
# # 	# limits = (0, nothing, 1, nothing)
# # 	# yticks = LogTicks(WilkinsonTicks(2))
# # 	# yticks = [10 ^ i for i in 1:9]
# # )
# # barplot!(ax, sorted_unique_distances, sorted_unique_distances_counts_cumsum, color = :red, strokecolor = :black, strokewidth = 1)
# # f

# f = Figure(size = (500, 500))
# ax = Axis(f[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Num of unique distances to chosen proteins per unchosen protein",
# 	ylabel = "Unchosen proteins (cumulative) [%]",
# 	# yscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	limits = (0, nothing, 0, nothing),
# 	xticks = 0:5:sorted_unique_distances[end],
# 	yticks = 0:10:100,
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# barplot!(ax, sorted_unique_distances, sorted_unique_distances_counts_cumsum_percent, color = :red, strokecolor = :black, strokewidth = 1)
# f



# # 2 - avg distance to chosen proteins

# avg_distance_per_unchosen_prot = mean.(eachrow(unchosenprot_chosenprot_distances))
# describe(avg_distance_per_unchosen_prot)


# f = Figure(size = (500, 500))
# ax = Axis(f[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Avg distance to chosen proteins per unchosen protein",
# 	ylabel = "Unchosen proteins",
# 	yscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# hist!(ax, avg_distance_per_unchosen_prot, color = :red, strokecolor = :black, strokewidth = 1)
# f


# # 3 - the 3 smallest distances between each unchosen protein to the chosen proteins
# # 4 - the number of chosen proteins per the 3 smallest distances for each unchosen protein


# function retain_smallest_distances_counter(unchosen_distances_counter::Accumulator, x_smallest::Int = 3)
# 	sorted_unchosen_distances = sort(collect(keys(unchosen_distances_counter)))
# 	smallest_distances = sorted_unchosen_distances[1:min(x_smallest, length(sorted_unchosen_distances))]
# 	# smallest_distances_counter = counter(Dict(d => unchosen_distances_counter[d] for d in smallest_distances))
# 	smallest_distances_counter = OrderedDict(d => unchosen_distances_counter[d] for d in smallest_distances)
# 	return smallest_distances_counter
# end

# smallest_unchosen_distances_counters = retain_smallest_distances_counter.(unchosen_distances_counters)

# # smallest_unchosen_distances_counters[1]

# f = Figure(size = (800, 400))
# axes = []
# for i in 1:3
# 	ax = Axis(f[1, i],
# 		yscale = log10,
# 		title = "Distance rank: $i",
# 		# limits = (0, nothing, 0, nothing),
# 		# limits = (0, nothing, 1, nothing),
# 		# yticks = range(0, 1; step=0.1),
# 	)
# 	push!(axes, ax)
# 	xs = []
# 	for smallest_unchosen_distances_counter in smallest_unchosen_distances_counters
# 		length(smallest_unchosen_distances_counter) < i && continue
# 		d = collect(keys(smallest_unchosen_distances_counter))[i]
# 		push!(xs, d)
# 	end
# 	hist!(
# 		ax,
# 		xs,
# 		color = :red,
# 		strokecolor = :black,
# 		strokewidth = 1,
# 		#  normalization = :probability
# 	)
# end
# linkxaxes!(axes...)
# linkyaxes!(axes...)
# Label(f[2, begin:end], "Distance to chosen protein")  # x axis title
# Label(f[begin:end, 0], "Unchosen proteins", rotation = pi / 2)  # y axis title
# Label(f[0, begin:end], "The 3 smallest distances between each unchosen protein to the chosen proteins", fontsize = 20)  # main title
# f


# f = Figure(size = (800, 400))
# all_xs = []
# for i in 1:3
# 	xs = []
# 	for smallest_unchosen_distances_counter in smallest_unchosen_distances_counters
# 		length(smallest_unchosen_distances_counter) < i && continue
# 		d = collect(keys(smallest_unchosen_distances_counter))[i]
# 		c = smallest_unchosen_distances_counter[d]
# 		push!(xs, c)
# 	end
# 	push!(all_xs, xs)
# end
# # max_x = maximum(maximum.(all_xs))
# axes = []
# for (i, xs) in enumerate(all_xs)
# 	ax = Axis(f[1, i],
# 		yscale = log10,
# 		title = "Distance rank: $i",
# 		# limits = (0, nothing, 0, nothing),
# 		limits = (0, nothing, 1, nothing),
# 		# xticks = 0:2500:max_x,
# 		# yticks = range(0, 1; step=0.1),
# 	)
# 	push!(axes, ax)
# 	hist!(
# 		ax,
# 		xs,
# 		color = :red,
# 		strokecolor = :black,
# 		strokewidth = 1,
# 		#  normalization = :probability
# 	)
# end
# linkxaxes!(axes...)
# linkyaxes!(axes...)
# Label(f[2, begin:end], "Chosen proteins in distance rank")  # x axis title
# Label(f[begin:end, 0], "Unchosen proteins", rotation = pi / 2)  # y axis title
# Label(f[0, begin:end], "The number of chosen prots per the 3 smallest distances for each unchosen prot", fontsize = 20)  # main title
# f




# principle_nas_per_unchosen = unchosendf[!, "MeanNAPositions"]
# all_xs = []
# all_nas_per_unchosen = []
# dfs = []
# for i in 1:3
# 	xs = []
# 	nas_per_unchosen = []
# 	for (smallest_unchosen_distances_counter, na_positions) in zip(smallest_unchosen_distances_counters, principle_nas_per_unchosen)
# 		length(smallest_unchosen_distances_counter) < i && continue
# 		d = collect(keys(smallest_unchosen_distances_counter))[i]
# 		c = smallest_unchosen_distances_counter[d]
# 		push!(xs, c)
# 		push!(nas_per_unchosen, na_positions)
# 	end
# 	push!(all_xs, xs)
# 	push!(all_nas_per_unchosen, nas_per_unchosen)
# 	df = DataFrame(
# 		"DistanceRank" => fill(i, length(xs)),
# 		"ChosenProteins" => xs,
# 		"MeanNAPositions" => nas_per_unchosen,
# 	)
# 	push!(dfs, df)
# end
# # max_x = maximum(maximum.(all_xs))

# f = Figure(size = (800, 400))
# axes = []
# for (i, (xs, nas_per_unchosen)) in enumerate(zip(all_xs, all_nas_per_unchosen))
# 	ax = Axis(f[1, i],
# 		# xscale = log10,
# 		yscale = log10,
# 		title = "Distance rank: $i",
# 		# limits = (0, nothing, 0, nothing),
# 		# limits = (0, nothing, 1, nothing),
# 		xticks = 0:20:maximum(principle_nas_per_unchosen),
# 		# yticks = range(0, 1; step=0.1),
# 	)
# 	push!(axes, ax)
# 	scatter!(
# 		ax,
# 		nas_per_unchosen,
# 		xs,
# 		# color = :red, strokecolor = :black, strokewidth = 1,
# 		#  normalization = :probability
# 	)
# end
# linkxaxes!(axes...)
# linkyaxes!(axes...)
# Label(f[2, begin:end], "Mean NA positions / unchosen protein")  # x axis title
# Label(f[begin:end-1, 0], "Chosen proteins", rotation = pi / 2)  # y axis title
# Label(f[0, begin+1:end], "Out-reassignments per unchosen protein vs. its reads' NAs", fontsize = 20)  # main title
# rowsize!(f.layout, 1, Relative(3 / 4))
# f



# avg_dfs = [
# 	combine(
# 		groupby(df, "MeanNAPositions"),
# 		nrow,
# 		:ChosenProteins => mean => :mean,
# 		:ChosenProteins => std => :std,
# 	)
# 	for df in dfs
# ]
# for avg_df in avg_dfs
# 	avg_df[!, "std"] = ifelse.(isnan.(avg_df[!, "std"]), 0, avg_df[!, "std"])
# end



# f = Figure(size = (800, 400))
# xtick_0 = 0
# xtick_d = 20
# xtick_upper_bound = maximum(principle_nas_per_unchosen) + xtick_d
# xticks = [
# 	xtick_0 + (xtick_d * i)
# 	for i in 0:div(xtick_upper_bound - xtick_0, xtick_d)+1
# ]
# if xticks[end] == xtick_upper_bound
# 	xticks = xticks[begin:end-1]
# end
# max_y = maximum([maximum(avg_df[!, "mean"] .+ avg_df[!, "std"]) for avg_df in avg_dfs])
# max_x = round(Int, max_y)
# max_y_magnitude = length(digits(max_x))
# yticks = 0:2000:10^max_y_magnitude
# axes = []
# for (i, df) in enumerate(avg_dfs)
# 	xs = df[!, "MeanNAPositions"]
# 	ys = df[!, "mean"]
# 	yerrs = df[!, "std"]
# 	# where y - yerr < 0, we want to make yerr = y
# 	# higherrors = copy(yerrs)
# 	# lowerrors = copy(yerrs)
# 	# lowerrors[ys.-lowerrors.<0] .= ys[ys.-lowerrors.<0]
# 	ax = Axis(f[1, i],
# 		# xscale = log10,
# 		# yscale = log10,
# 		title = "Distance rank: $i",
# 		# limits = (0, nothing, 0, nothing),
# 		# limits = (0, nothing, 1, nothing),
# 		xticks = xticks,
# 		yticks = yticks,
# 		# yminorgridvisible = true,
# 		# yminorticks = IntervalsBetween(5)
# 	)
# 	# ax.yticks = 0:2000:
# 	push!(axes, ax)
# 	errorbars!(
# 		ax,
# 		xs,
# 		ys,
# 		yerrs,
# 		# lowerrors, higherrors,
# 		color = :red,
# 		label = :std,
# 	)
# 	# plot position scatters so low and high errors can be discriminated
# 	scatter!(ax, xs, ys; markersize = 4, color = :black, label = :mean)
# end
# linkxaxes!(axes...)
# linkyaxes!(axes...)
# Label(f[2, begin:end], "Mean NA positions / unchosen protein")  # x axis title
# Label(f[begin:end-1, 0], "Mean chosen proteins", rotation = pi / 2)  # y axis title
# Label(f[0, begin+1:end], "Out-reassignments per unchosen protein vs. its reads' NAs", fontsize = 20)  # main title
# rowsize!(f.layout, 1, Relative(3 / 4))
# f



# c_or_unc_titles = ["Chosen", "Unchosen"]
# f = Figure(size = (800, 400))
# axes = []
# for (i, (df, title)) in enumerate(zip(c_or_unc_dfs, c_or_unc_titles))
# 	xs = df[!, "MeanNAPositions"]
# 	ax = Axis(f[1, i],
# 		# xscale = log10,
# 		# yscale = log10,
# 		title = title,
# 		# limits = (0, nothing, 0, nothing),
# 		# limits = (0, nothing, 1, nothing),
# 		# xticks = xticks,
# 		# yticks = yticks,
# 		# yminorgridvisible = true,
# 		# yminorticks = IntervalsBetween(5)
# 	)
# 	# ax.yticks = 0:2000:
# 	push!(axes, ax)
# 	hist!(
# 		ax,
# 		xs,
# 		color = :red,
# 		strokecolor = :black,
# 		strokewidth = 1,
# 		normalization = :probability,
# 	)
# end
# linkxaxes!(axes...)
# linkyaxes!(axes...)
# Label(f[2, begin:end], "Mean NA positions")  # x axis title
# # Label(f[begin:end, 0], "Proteins", rotation = pi / 2)  # y axis title
# Label(f[begin:end-1, 0], "Proteins (probability)", rotation = pi / 2, tellheight = false)  # y axis title
# # Label(f[0, begin+1:end], "Reassignments vs. NAs", fontsize = 20)  # main title
# f








# one_result = one_solution_additional_assignment_considering_available_reads(
# 	distinctdf, allprotsdf, firstcolpos, Δ, solution,
# 	readsdf
# )






# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	xlabel = "Mean NA positions / chosen protein",
# 	ylabel = "Mean reassigned unchosen proteins",
# 	# yscale = log10,
# 	title = "In-reassignments per chosen protein vs. its reads' NAs",
# 	# limits = (0, nothing, 1, nothing)
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# df = combine(
# 	groupby(one_result, "MeanNAPositions"),
# 	nrow,
# 	:AdditionalSupportingProteins => mean => :mean,
# 	:AdditionalSupportingProteins => std => :std,
# )
# xs = df[!, "MeanNAPositions"]
# ys = df[!, "mean"]
# yerrs = df[!, "std"]
# errorbars!(
# 	ax,
# 	xs,
# 	ys,
# 	yerrs,
# 	# lowerrors, higherrors,
# 	color = :red,
# 	# label = :std,
# )
# # plot position scatters so low and high errors can be discriminated
# scatter!(ax, xs, ys; markersize = 8, color = :black)
# fig




# # df = deepcopy(one_result)
# # df[!, "MeanNAPositionsPerAdditionalSupportingProtein"] = mean.(one_result[!, "AdditionalSupportingProteinsMeanNAPositions"])
# # df = combine(
# # 	groupby(df, "MeanNAPositions"),
# # 	nrow,
# # 	:MeanNAPositionsPerAdditionalSupportingProtein => mean => :mean,
# # 	:MeanNAPositionsPerAdditionalSupportingProtein => std => :std,
# # )
# # fig = Figure(size = (500, 500))
# # ax = Axis(fig[1, 1],
# # 	xlabel = "Mean NA positions / chosen protein",
# # 	ylabel = "Mean mean NA positions / reassigned unchosen protein",
# # 	# yscale = log10,
# # 	# title = "In-reassignments per chosen protein vs. its reads' NAs",
# # 	# limits = (0, nothing, 1, nothing)
# # 	# yticks = LogTicks(WilkinsonTicks(2))
# # 	# yticks = [10 ^ i for i in 1:9]
# # )
# # xs = df[!, "MeanNAPositions"]
# # ys = df[!, "mean"]
# # yerrs = df[!, "std"]
# # errorbars!(
# # 	ax,
# # 	xs,
# # 	ys,
# # 	yerrs,
# # 	# lowerrors, higherrors,
# # 	color = :red,
# # 	# label = :std,
# # )
# # # plot position scatters so low and high errors can be discriminated
# # scatter!(ax, xs, ys; markersize = 8, color = :black)
# # fig




# df = flatten(
# 	one_result[!, ["MeanNAPositions", "AdditionalSupportingProteinsMeanNAPositions"]],
# 	"AdditionalSupportingProteinsMeanNAPositions"
# )
# rename!(df, "AdditionalSupportingProteinsMeanNAPositions" => "AdditionalSupportingProteinMeanNAPositions")
# df = combine(
# 	groupby(df, "MeanNAPositions"),
# 	nrow,
# 	:AdditionalSupportingProteinMeanNAPositions => mean => :mean,
# 	:AdditionalSupportingProteinMeanNAPositions => std => :std,
# )
# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	xlabel = "Mean NA positions / chosen protein",
# 	ylabel = "Mean NA positions / reassigned unchosen protein",
# 	# yscale = log10,
# 	# title = "In-reassignments per chosen protein vs. its reads' NAs",
# 	# limits = (0, nothing, 1, nothing)
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# xs = df[!, "MeanNAPositions"]
# ys = df[!, "mean"]
# yerrs = df[!, "std"]
# errorbars!(
# 	ax,
# 	xs,
# 	ys,
# 	yerrs,
# 	# lowerrors, higherrors,
# 	color = :red,
# 	# label = :std,
# )
# # plot position scatters so low and high errors can be discriminated
# scatter!(ax, xs, ys; markersize = 8, color = :black)
# fig





# sum(length.(split.(one_result[!, "AdditionalSupportingProteinsIDs"], ",")))


# reassigned_prots_counter = counter(vcat(split.(one_result[!, "AdditionalSupportingProteinsIDs"], ",")...))
# reassignments_per_reassigned_prot = collect(values(reassigned_prots_counter))
# describe(reassignments_per_reassigned_prot)

# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Reassignments to X chosen proteins",
# 	ylabel = "Unchosen proteins",
# 	yscale = log10,
# 	xscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]
# )
# hist!(
# 	ax, reassignments_per_reassigned_prot, color = :red, 
# 	strokecolor = :red, strokewidth = 1,
# 	bins=maximum(reassignments_per_reassigned_prot),
# )
# fig



# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Reassignments to X chosen proteins",
# 	ylabel = "Unchosen proteins (cumulative prob)",
# 	xscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (1, nothing, 0, 1),
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	yticks = 0:0.1:1
# )
# ecdfplot!(
# 	ax, reassignments_per_reassigned_prot, 
# 	# color = :red, 
# 	# strokecolor = :black, strokewidth = 1,
# 	# bins=maximum(reassignments_per_reassigned_prot),
# )
# fig



# df = deepcopy(one_result)
# df[!, "%OriginalReads/TotalReads"] = 100 .* df[!, "NumOfReads"] ./ df[!, "TotalWeightedSupportingReads"]
# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	# xlabel = "Distance", 

# 	xlabel = "Original reads / total reads [%]",
# 	ylabel = "Chosen proteins [probability]",
# 	# yscale = log10,
# 	# xscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	xticks=0:10:100,
# 	# yticks = LogTicks(WilkinsonTicks(2))
# 	# yticks = [10 ^ i for i in 1:9]

# )
# hist!(
# 	ax, df[!, "%OriginalReads/TotalReads"], color = :red, 
# 	strokecolor = :black, strokewidth = 1,
# 	bins=20,
# 	normalization = :probability
# )
# fig


# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	xlabel = "Original reads / chosen protein",
# 	ylabel = "Reassigned unchosen protein",
# 	xscale = log10,
# 	yscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	# xticks=0:10:100,
# 	# yticks = [10 ^ i for i in 1:9]
# )
# xs = one_result[!, "NumOfReads"]
# ys = one_result[!, "AdditionalSupportingProteins"]
# scatter!(
# 	ax, xs, ys 
# 	# color = :red, 
# 	# strokecolor = :black, strokewidth = 1,
# 	# bins=20,
# 	# normalization = :probability
# )
# fig


# fig = Figure(size = (500, 500))
# ax = Axis(fig[1, 1],
# 	xlabel = "Original reads / chosen protein",
# 	ylabel = "Total reads",
# 	xscale = log10,
# 	yscale = log10,
# 	# title = "Reassignments per unchosen protein",
# 	# limits = (0, nothing, 1, nothing)
# 	# xticks=0:10:100,
# 	# yticks = [10 ^ i for i in 1:9]
# )
# xs = one_result[!, "NumOfReads"]
# ys = one_result[!, "TotalWeightedSupportingReads"]
# scatter!(
# 	ax, xs, ys 
# 	# color = :red, 
# 	# strokecolor = :black, strokewidth = 1,
# 	# bins=20,
# 	# normalization = :probability
# )
# fig

































































# unique(finalresults[!, "#Solution"])


# distinctfiles = split(readchomp(`cat "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/DistinctProteinsForExpressionLevels.txt"`))
# allprotsfiles = split(readchomp(`cat "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/UniqueProteinsForExpressionLevels.txt"`))
# samplenames = split(readchomp(`cat "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/ChromsNamesForExpressionLevels.txt"`))

# `echo O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/*.ExpressionLevels.csv`

# run(`echo "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/*.ExpressionLevels.csv"`)

# read(`echo "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/*.ExpressionLevels.csv"`, String)

# expressionfiles = split(read(`echo "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/*.ExpressionLevels.csv"`, String))

# protsdir = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles"

# expressionfiles = [
#     f for f in readdir(protsdir) if occursin("ExpressionLevels.csv", f)
# ]

# outdir = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles"
# postfix_to_add = ""
# firstcolpos = 16
# delim = "\t"
# innerdelim = ","
# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]
# # fractions = [1.0]
# fractions = [0.2, 0.4, 0.6, 0.8, 1.0]
# maxmainthreads = 30
# algs = ["Ascending", "Descending"]
# onlymaxdistinct = false
# gcp = false
# shutdowngcp = false
# substitutionmatrix = nothing
# aagroups = nothing
# similarityscorecutoff = 0
# similarityvalidator = >=

# x = 59

# distinctdf = prepare_distinctdf(
#     distinctfiles[x], delim, innerdelim, truestrings, falsestrings
# )

# allprotsdf, firstcolpos = prepare_allprotsdf!(
#     allprotsfiles[x], delim, innerdelim, truestrings, falsestrings, firstcolpos
# )

# M = Matrix(allprotsdf[:, firstcolpos:end])

# Δ = distances(M)

# solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)


# main(
#     distinctfiles, allprotsfiles, samplenames, postfix_to_add,
#     firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
#     maxmainthreads, outdir, algs, onlymaxdistinct,
#     gcp, shutdowngcp,
#     substitutionmatrix,
#     similarityscorecutoff,
#     similarityvalidator,
#     aagroups
# )





# distinctfiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv"
# ]
# allprotsfiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# ]
# samplenames = [
#     "GRIA",
# ]
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3"


# # input

# distinctfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.03.03.2023-15:36:38.csv"
# allprotsfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# samplename = "GRIA"
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3"
# postfix_to_add = ""

# firstcolpos = 15
# delim = "\t"
# innerdelim = ","
# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]
# fractions = [1.0]
# maxmainthreads = 30
# algs = ["Ascending", "Descending"]
# onlymaxdistinct = false


# # # # run_sample

# distinctdf = prepare_distinctdf(
#     distinctfile, delim, innerdelim, truestrings, falsestrings
# )

# # considering only desired solutions (rows' indices)
# # solutions = distinctdf[:, "Index"]

# _distinctdf1 = subset(
#     distinctdf, 
#     "Fraction" => x -> x .∈ fractions,  # keep only solutions of desired fractions
#     # "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
# )

# occursin.("Ascending", Ref(algs))
# occursin.("Ascending", algs)
# any(occursin.("Ascending", algs))

# _distinctdf2 = subset(
#     distinctdf, 
#     # "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
#     "Algorithm" => x -> occursin.(x, Ref(algs)) # keep only solutions of desired algortihms
#     # "Algorithm" => x -> any(occursin.(x, algs)) # keep only solutions of desired algortihms
# )

# subset(
#     distinctdf[!, ["Algorithm"]],
#     # "Algorithm" => x -> any(occursin.(x, algs)) # keep only solutions of desired algortihms
#     "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
#     # "Algorithm" => x -> x .∈ algs # keep only solutions of desired algortihms
# )




# df = DataFrame(
#     Algorithm = ["Ascending", "Ascending", "Naive", "Naive", "Descending", "Descending"],
#     Fraction = [0.5, 1.0, 0.5, 1.0, 0.5, 1.0]
# )
# algs = ["Ascending", "Descending"]
# fractions = [1.0]
# # subset(
# #     df,
# #     "Algorithm" => x -> any(occursin.(x, algs)) # keep only solutions of desired algortihms
# #     # "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
# #     # "Algorithm" => x -> x .∈ algs # keep only solutions of desired algortihms
# # )

# in.(df.Algorithm, Ref(algs))
# in.(df.Fraction, Ref(fractions))

# in.(df.Algorithm, Ref(algs)) .& in.(df.Fraction, Ref(fractions))

# df[in.(df.Algorithm, Ref(algs)) .& in.(df.Fraction, Ref(fractions)), :]

# occursin.("Ascending", algs)


# a = [true, false]
# b = [true, true]
# a .& b
# a .| b

# _distinctdf = subset(
#     distinctdf, 
#     "Fraction" => x -> x .∈ fractions,  # keep only solutions of desired fractions
#     "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
# )


# solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)

# # considering only solutions of desired fractions

# algs = ["Ascending"]


# distinctdf = subset(distinctdf, "Fraction" => x -> x .∈ fractions, "Algorithm" => x -> occursin.(x, algs))


# maxdistinct = maximum(distinctdf[!, "NumUniqueSamples"])
# distinctdf = subset(distinctdf, "NumUniqueSamples" => x -> x .== maxdistinct)



# subset(distinctdf, "Fraction" => x -> x .∈ fractions, "NumUniqueSamples" => x -> x .== maxdistinct)


# # original
# solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[:, "Index"]


# # basesize = Int(round(Threads.nthreads()))
# # maxmainthreads = Int(round(Threads.nthreads() / 5))
# # basesize = Int(round(Threads.nthreads() / 5))
# allsubsolutions = collect(Iterators.partition(solutions, maxmainthreads))
# results = tcollect(
#     additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
#     for subsolutions ∈ allsubsolutions
# )
# finalresults = vcat(Iterators.flatten(results)...)

# # save the results
# outfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels$postfix_to_add.csv")
# CSV.write(outfile, finalresults; delim)















# main(
#     distinctfiles, allprotsfiles, samplenames,
#     firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
#     maxmainthreads, outdir
# )


