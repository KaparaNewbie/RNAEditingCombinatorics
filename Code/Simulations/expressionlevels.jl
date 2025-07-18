using ArgParse
using CSV
using DataFrames
using StatsBase  # for StatsBase.sample
using BenchmarkTools
using BioSequences # for BioSequences.toAAset
import Base.Threads.@spawn
using ThreadsX # for ThreadsX.Set, 
using Transducers  # for tcollect
using SymmetricFormats # for SymmetricPacked which creates a symmetric matrix from a triangular array
using BioAlignments # for SubstitutionMatrix


include(joinpath(@__DIR__, "consts.jl")) # for `AA_groups` (polar/non-polar/positive/negative)
include(joinpath(@__DIR__, "issimilar.jl")) # for issimilar
include(joinpath(@__DIR__, "timeformatters.jl"))



"""
	distances(M, aagroups)

Create a symmaetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two contain no possible amino acids `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, 
such that both `AAᵦ` and `AAᵧ` share the same classification in `aagroups`.
"""
function distances(M::Matrix{Set{AminoAcid}}, aagroups::Dict{AminoAcid, String})
	nrows = size(M, 1)
	Δmax = size(M, 2)  # maximal possible difference between any two proteins
	uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	# sets with no two possible AAs that share a similar classification
	distinctAAsets = ThreadsX.Set(
		[
		(x, y)
		for x ∈ AAsets for y ∈ AAsets
		if !anysimilarity(x, y, aagroups)
	]
	)
	# for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
	Δs = tcollect(i_js_distances(M, distinctAAsets, uint, i) for i ∈ 1:nrows)
	# merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
	Δ = SymmetricPacked(reduce(vcat, Δs))
	return Δ
end




"""
	distances(M, substitutionmatrix, similarityscorecutoff, similarityvalidator)

Create a symmaetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acis `(AAᵦ, AAᵧ) ∈ (Sᵢ x Sⱼ)`, 
such that their substitution score according to `substitutionmatrix` is `>`/`≥`/`<`/`≤` `similarityscorecutoff`. 
The specific comparison (e.g., `≥`) is determined by `similarityvalidator`.
"""
function distances(
	M::Matrix{Set{AminoAcid}},
	substitutionmatrix::SubstitutionMatrix{AminoAcid, Int64}, similarityscorecutoff::Int64, similarityvalidator::Function,
)
	nrows = size(M, 1)
	Δmax = size(M, 2)  # maximal possible difference between any two proteins
	uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	# sets with no two possible AAs whose substitution score is acceptable
	distinctAAsets = ThreadsX.Set(
		[
		(x, y)
		for x ∈ AAsets for y ∈ AAsets
		if !anysimilarity(x, y, substitutionmatrix, similarityscorecutoff, similarityvalidator)
	]
	)
	# for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
	Δs = tcollect(i_js_distances(M, distinctAAsets, uint, i) for i ∈ 1:nrows)
	# merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
	Δ = SymmetricPacked(reduce(vcat, Δs))
	return Δ
end




"""
	distances(M)

Create a symmaetrical distances matrix `Δ` which measures the distance between any `rowᵢ, rowⱼ ∈ M`.  
The distance between a `rowᵢ` to a `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acid.
"""
function distances(M::Matrix{Set{AminoAcid}})
	nrows = size(M, 1)
	Δmax = size(M, 2)  # maximal possible difference between any two proteins
	uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	# sets whose intersection is empty
	distinctAAsets = ThreadsX.Set(
		[(x, y)
		 for x ∈ AAsets
		 for y ∈ AAsets
		 if x ∩ y == Set()]
	)
	# for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
	Δs = tcollect(i_js_distances(M, distinctAAsets, uint, i) for i ∈ 1:nrows)
	# merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
	Δ = SymmetricPacked(reduce(vcat, Δs))
	return Δ
end


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
	@views rowᵢ = M[i, :]
	nrows = size(M, 1)
	δᵢ = Vector{uint}(undef, nrows - i + 1)
	# the distance between `rowᵢ` to itself
	δᵢ[1] = 0
	i == nrows && return δᵢ
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
	return δᵢ
end



const UInts = [UInt8, UInt16, UInt32, UInt64, UInt128]


function smallestuint(maximumval)
	for uint ∈ UInts
		if maximumval <= typemax(uint)
			return uint
		end
		return UInts[end]
	end
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




function run_sample(
	distinctfile, allprotsfile, samplename, postfix_to_add,
	firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
	maxmainthreads, outdir, algs, onlymaxdistinct,
	substitutionmatrix::Union{SubstitutionMatrix, Nothing},
	similarityscorecutoff::Int64,
	similarityvalidator::Function,
	aagroups::Union{Dict{AminoAcid, String}, Nothing},
)
	@info "$(loggingtime())\trun_sample" distinctfile allprotsfile samplename

	distinctdf = prepare_distinctdf(
		distinctfile, delim, innerdelim, truestrings, falsestrings,
	)

	allprotsdf, firstcolpos = prepare_allprotsdf!(
		allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos,
	)

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

	minmainthreads = minimum([Int(Threads.nthreads() / 5), Int(length(solutions) / 4)])
	allsubsolutions = collect(Iterators.partition(solutions, minmainthreads))

	results = tcollect(
		additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
		for subsolutions ∈ allsubsolutions
	)
	# finalresults = vcat(Iterators.flatten(results)...)
	# finalresults = vcat(skipmissing(Iterators.flatten(results)...))
	finalresults = vcat((skipmissing(Iterators.flatten(results))...))

	# save the results
	outfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels$postfix_to_add.csv")
	CSV.write(outfile, finalresults; delim)
end


"Get indices of eligible solutions from `distinctdf` on which expression levels should be calculated."
function choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)
	# _distinctdf = subset(
	#     distinctdf, 
	#     "Fraction" => x -> x .∈ fractions,  # keep only solutions of desired fractions
	#     "Algorithm" => x -> occursin.(x, algs) # keep only solutions of desired algortihms
	# )
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




function additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, solutions)
	results = map(solutions) do solution
		try
			one_solution_additional_assignment_considering_available_reads(distinctdf, allprotsdf, firstcolpos, Δ, solution)
		catch e
			@warn "Error in additional_assignments:" e solution
			missing
		end
	end
	return results
end



function prepare_distinctdf(
	distinctfile, delim, innerdelim,
	truestrings,
	# truestrings::Union{Nothing, Vector{String}, Vector{Any}},
	falsestrings,
	# falsestrings::Union{Nothing, Vector{String}, Vector{Any}},
)
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


# firstcolpos = 15
function prepare_allprotsdf!(allprotsfile, delim, innerdelim,
	truestrings,
	# truestrings::Union{Nothing, Vector{String}, Vector{any}},
	falsestrings, firstcolpos,
)
	# allprotsfile = allprotsfiles[1]

	# allprotsdf = DataFrame(CSV.File(allprotsfile; delim, truestrings, falsestrings))
	# allprotsdf = hcat(allprotsdf[:, begin:firstcolpos-1], toAAset.(allprotsdf[:, firstcolpos:end], innerdelim))

	df1 = DataFrame(CSV.File(allprotsfile, delim = delim, select = collect(1:firstcolpos-1), types = Dict("Protein" => String, "Reads" => String)))
	df1[!, "Protein"] = InlineString.(df1[!, :Protein])
	# make sure columns of AAs containing only Ts aren't parsed as boolean columns
	df2 = DataFrame(CSV.File(allprotsfile, delim = delim, drop = collect(1:firstcolpos-1), types = String))
	df2 = toAAset.(df2, innerdelim)

	allprotsdf = hcat(df1, df2)

	transform!(allprotsdf, :Reads => (x -> split.(x, innerdelim)) => :Reads)

	# typeof(allprotsdf[!, :Reads][1][1])
	# vec{SubString{String}}()
	# eltype
	# vec{typeof(allprotsdf[!, :Reads][1])}()

	# readtype = typeof(allprotsdf[!, :Reads][1][1])
	# Vector{readtype}(undef, 0)
	# Matrix{readtype}(undef, size(allprotsdf, 1), 0)
	# [Vector{readtype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)]

	insertcols!(allprotsdf, firstcolpos, :Index => 1:size(allprotsdf, 1))
	firstcolpos += 1
	insertcols!(allprotsdf, firstcolpos, :AdditionalEqualSupportingReads => 0.0)
	firstcolpos += 1
	insertcols!(allprotsdf, firstcolpos, :AdditionalWeightedSupportingReads => 0.0)
	firstcolpos += 1

	# typeof(allprotsdf[!, :Reads][1])
	# eltype(allprotsdf[!, :Reads]) == typeof(allprotsdf[!, :Reads][1])
	# typeof(allprotsdf[!, :Reads][1][1])
	# readtype = typeof(allprotsdf[!, :Reads][1][1])
	# insertcols!(allprotsdf, firstcolpos, :AdditionalEqualSupportingReadsIDs => [Vector{readtype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1
	# insertcols!(allprotsdf, firstcolpos, :AdditionalWeightedSupportingReadsIDs => [Vector{readtype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1

	# readsvectype = eltype(allprotsdf[!, :Reads])
	# # insertcols!(allprotsdf, firstcolpos, :AdditionalEqualSupportingReadsIDs => [Vector{readsvectype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# # firstcolpos += 1
	# # insertcols!(allprotsdf, firstcolpos, :AdditionalWeightedSupportingReadsIDs => [Vector{readsvectype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# # firstcolpos += 1
	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingReadsIDs => [Vector{readsvectype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1
	# proteinseltype = eltype(allprotsdf[!, :Protein])
	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingProteinsIDs => [Vector{proteinseltype}(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1


	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingReadsIDs => [Vector(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1
	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingProteinsIDs => [Vector(undef, 0) for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1

	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingReadsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1
	# insertcols!(allprotsdf, firstcolpos, :AdditionalSupportingProteinsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)])
	# firstcolpos += 1

	insertcols!(
		allprotsdf,
		firstcolpos,
		:AdditionalSupportingReadsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)],
		:AdditionalSupportingProteinsIDs => [[] for _ ∈ 1:size(allprotsdf, 1)],
	)
	firstcolpos += 2

	# allprotsdf[:, firstcolpos-6:firstcolpos]

	# transform!(allprotsdf, :Reads => (x -> split.(x, innerdelim)) => :Reads)

	return allprotsdf, firstcolpos
end


# solution = 901

# unchosendf = one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, solution)

# describe(unchosendf[!, ["NumOfReads", "AdditionalEqualSupportingReads", "AdditionalWeightedSupportingReads", "TotalEqualSupportingRead", "TotalWeightedSupportingRead"]])
# sum.(eachcol(unchosendf[!, ["NumOfReads", "AdditionalEqualSupportingReads", "AdditionalWeightedSupportingReads", "TotalEqualSupportingRead", "TotalWeightedSupportingRead"]]))

# function one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, solution)

#     solutionrow = distinctdf[solution, :]

#     prots_in_solution = solutionrow["UniqueSamples"]

#     allprotsdf = allprotsdf[:, begin:firstcolpos-1]

#     # sum(allprotsdf[:, "NumOfReads"])

#     chosendf = filter("Protein" => x -> x ∈ prots_in_solution, allprotsdf)
#     unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, allprotsdf)

#     insertcols!(chosendf, 3, "#Solution" => solutionrow["Index"])
#     insertcols!(chosendf, 4, "Fraction" => solutionrow["Fraction"])
#     insertcols!(chosendf, 5, "FractionRepetition" => solutionrow["FractionRepetition"])
#     insertcols!(chosendf, 6, "Algorithm" => solutionrow["Algorithm"])
#     insertcols!(chosendf, 7, "AlgorithmRepetition" => solutionrow["AlgorithmRepetition"])

#     chosenindices = chosendf[!, "Index"] # indices of chosen proteins in the complete Δ matrix

#     for unchosenprot ∈ eachrow(unchosendf)  # a row of unchosen protein relative to the chosen proteins in the solution

#         # unchosenprot = unchosendf[1000, :]  # an example row of unchosen protein

#         unchosenprot_index = unchosenprot["Index"] # the row number of the unchosen protein in the distances matrix Δ

#         unchosenprot_distances = Δ[unchosenprot_index, chosenindices] # the distances between the unchosen protein `unchosenprot` to all chosen proteins

#         unchosenprot_distances_argmins = allargmins(unchosenprot_distances) # indices of minimum distances

#         minchosendf = chosendf[unchosenprot_distances_argmins, :] # chosen proteins with minimum distance to the unchosen protein `unchosenprot`

#         unchosenreads = unchosenprot["NumOfReads"] # the number of reads `unchosenprot` divides between the chosen proteins that are closest to it

#         equal_addition = unchosenreads / size(minchosendf, 1)
#         weighted_addition = unchosenreads .* minchosendf[!, "NumOfReads"] ./ sum(minchosendf[!, "NumOfReads"])

#         # sum(chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"])
#         # sum(chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"])

#         chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"] .+= equal_addition
#         chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"] .+= weighted_addition

#         # sum(chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"])
#         # sum(chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"])

#     end

#     chosendf[!, "TotalEqualSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalEqualSupportingReads"]
#     chosendf[!, "TotalWeightedSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalWeightedSupportingReads"]

#     return chosendf

# end
# 



# firstcolpos = 15
# allprotsdf, firstcolpos = prepare_allprotsdf!(
#     allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos
# )
# allprotsdf[:, firstcolpos-4:firstcolpos-1]
# @assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])
# @assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])

# one_solution_additional_assignment_considering_available_reads(distinctdf, allprotsdf, firstcolpos, Δ, 7)
# one_solution_additional_assignment_considering_available_reads(distinctdf, allprotsdf, firstcolpos, Δ, 10)


function one_solution_additional_assignment_considering_available_reads(distinctdf, allprotsdf, firstcolpos, Δ, solution)

	# solution = 7 # todo comment out
	# solution = 10 # todo comment out

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
		baseallprotsdf = filter("NumOfReads" => x -> x > 0, baseallprotsdf)
	end

	# ["a1", "a2", "a3",  "b", "c"], Protein = ["A", "A", "A", "B", "C"]

	# availablereads = ["a1", "a3",  "b",]
	# allreadsperprotein = [["a1", "a2", "a3"],  ["b"], ["c"]]
	# availablereadsperprotein = [
	#         [read for read ∈ reads if read ∈ availablereads]
	#         for reads ∈ allreadsperprotein
	#     ]

	prots_in_solution = solutionrow["UniqueSamples"]

	chosendf = filter("Protein" => x -> x ∈ prots_in_solution, baseallprotsdf)
	unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, baseallprotsdf)

	insertcols!(
		chosendf,
		3,
		"#Solution" => solutionrow["Index"],
		"Fraction" => solutionrow["Fraction"],
		"FractionRepetition" => solutionrow["FractionRepetition"],
		"Algorithm" => solutionrow["Algorithm"],
		"AlgorithmRepetition" => solutionrow["AlgorithmRepetition"],
	)

	chosenindices = chosendf[:, "Index"] # indices of chosen proteins in the complete Δ matrix

	for unchosenprot ∈ eachrow(unchosendf)  # a row of unchosen protein relative to the chosen proteins in the solution

		# unchosenprot = unchosendf[1, :]  # TODO comment out: an example row of unchosen protein
		# unchosenprot = unchosendf[2, :]  

		unchosenprot_index = unchosenprot["Index"] # the row number of the unchosen protein in the distances matrix Δ

		unchosenprot_distances = Δ[unchosenprot_index, chosenindices] # the distances between the unchosen protein `unchosenprot` to all chosen proteins

		unchosenprot_distances_argmins = allargmins(unchosenprot_distances) # indices of minimum distances

		minchosendf = chosendf[unchosenprot_distances_argmins, :] # chosen proteins with minimum distance to the unchosen protein `unchosenprot`

		unchosenreads = unchosenprot["NumOfReads"] # the number of reads `unchosenprot` divides between the chosen proteins that are closest to it

		equal_addition = unchosenreads / size(minchosendf, 1)
		weighted_addition = unchosenreads .* minchosendf[:, "NumOfReads"] ./ sum(minchosendf[:, "NumOfReads"])

		chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"] .+= equal_addition
		chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"] .+= weighted_addition

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

	end

	# chosendf[!, "TotalEqualSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalEqualSupportingReads"]
	# chosendf[!, "TotalWeightedSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalWeightedSupportingReads"]
	# chosendf[!, "AdditionalSupportingProteins"] .= length.(chosendf[!, "AdditionalSupportingProteinsIDs"])

	chosendf[:, "TotalEqualSupportingReads"] .= chosendf[:, "NumOfReads"] .+ chosendf[:, "AdditionalEqualSupportingReads"]
	chosendf[:, "TotalWeightedSupportingReads"] .= chosendf[:, "NumOfReads"] .+ chosendf[:, "AdditionalWeightedSupportingReads"]
	chosendf[:, "AdditionalSupportingProteins"] .= length.(chosendf[:, "AdditionalSupportingProteinsIDs"])

	# @assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])
	# @assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])

	sort!(chosendf, "AdditionalSupportingProteins")

	# chosendf[:, 23:end]

	transform!(chosendf, :AdditionalSupportingProteinsIDs => ByRow(x -> join(x, ",")) => :AdditionalSupportingProteinsIDs)
	transform!(chosendf, :AdditionalSupportingReadsIDs => ByRow(x -> join(join.(x, ","), ";")) => :AdditionalSupportingReadsIDs) # todo uncomment after ensuring no undef values snaeak in


	@assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution"""
	@assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)]) """Solution: $solution"""

	# chosendf[:, "AdditionalSupportingProteinsIDs"]
	# # typeof.(chosendf[:, "AdditionalSupportingProteinsIDs"])
	# # isassigned.(chosendf[:, "AdditionalSupportingProteinsIDs"])
	# # chosendf[!, !][map(!, isassigned.(chosendf[!, "AdditionalSupportingProteinsIDs"]))]
	# # 
	# chosendf[:, "AdditionalSupportingReadsIDs"]

	# badchosendf = filter("AdditionalSupportingProteinsIDs" => x -> !isassigned(x), chosendf)
	# badchosendf[:, 23:end]
	# transform!(badchosendf, :AdditionalSupportingProteinsIDs => ByRow(x -> join(x, ",")) => :AdditionalSupportingProteinsIDs)
	# transform!(badchosendf, :AdditionalSupportingReadsIDs => ByRow(x -> join(join.(x, ","), ";")) => :AdditionalSupportingReadsIDs) # todo uncomment after ensuring no undef values snaeak in
	# badchosendf[:, 23:end]

	# goodchosendf = filter("AdditionalSupportingProteinsIDs" => x -> isassigned(x), chosendf)
	# goodchosendf[:, 23:end]
	# transform!(goodchosendf, :AdditionalSupportingProteinsIDs => ByRow(x -> join(x, ",")) => :AdditionalSupportingProteinsIDs)
	# goodchosendf[:, 23:end]
	# transform!(goodchosendf, :AdditionalSupportingReadsIDs => ByRow(x -> join(join.(x, ","), ";")) => :AdditionalSupportingReadsIDs) # todo uncomment after ensuring no undef values snaeak in
	# goodchosendf[:, 23:end]

	return chosendf

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

	end
	return parse_args(s)
end



function main(
	distinctfiles, allprotsfiles, samplenames,
	postfix_to_add,
	firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
	maxmainthreads, outdir, algs, onlymaxdistinct,
	gcp, shutdowngcp,
	substitutionmatrix::Union{SubstitutionMatrix, Nothing},
	similarityscorecutoff::Int64,
	similarityvalidator::Function,
	aagroups::Union{Dict{AminoAcid, String}, Nothing},
)
	@info "$(loggingtime())\tmain" distinctfiles allprotsfiles samplenames postfix_to_add firstcolpos delim innerdelim truestrings falsestrings fractions maxmainthreads outdir algs onlymaxdistinct gcp shutdowngcp substitutionmatrix similarityscorecutoff similarityvalidator aagroups

	# run each sample using the `run_sample` function
	for (distinctfile, allprotsfile, samplename) ∈ zip(distinctfiles, allprotsfiles, samplenames)
		run_sample(
			distinctfile, allprotsfile, samplename, postfix_to_add,
			firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
			maxmainthreads, outdir, algs, onlymaxdistinct,
			substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups,
		)
	end

	# shutdown gcp vm (if needed)
	gcp && shutdowngcp && run(`sudo shutdown`) # https://cloud.google.com/compute/docs/shutdownscript
end


function CLI_main()
	# read command-line args
	parsedargs = parsecmd()

	distinctfiles = parsedargs["distinctfiles"]
	allprotsfiles = parsedargs["allprotsfiles"]
	samplenames = parsedargs["samplenames"]

	distinctfilesfofn = parsedargs["distinctfilesfofn"]
	allprotsfilesfofn = parsedargs["allprotsfilesfofn"]
	samplenamesfile = parsedargs["samplenamesfile"]

	postfix_to_add = parsedargs["postfix_to_add"]

	firstcolpos = parsedargs["firstcolpos"]
	delim = parsedargs["delim"]
	innerdelim = parsedargs["innerdelim"]
	truestrings = parsedargs["truestrings"]
	falsestrings = parsedargs["falsestrings"]

	fractions = parsedargs["fractions"]
	fractions = fractions isa Vector{Float64} ? fractions : parse.(Float64, fractions)
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

	if distinctfilesfofn != ""
		distinctfiles = split(readline(distinctfilesfofn), " ")
	end
	if allprotsfilesfofn != ""
		allprotsfiles = split(readline(allprotsfilesfofn), " ")
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
		substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups,
	)
end


if abspath(PROGRAM_FILE) == @__FILE__
	CLI_main()
end




# # distinctfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.DistinctUniqueProteins.Fraction0_1.11.05.2023-17:03:49.csv"
# distinctfile = "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/DistinctProteins/comp73852_c0_seq1.DistinctUniqueProteins.29.01.2025-14:49:19.csv"
# delim = "\t"
# innerdelim = ","
# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]
# # allprotsfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# allprotsfile = "/private7/projects/Combinatorics/O.vulgaris/MpileupAndTranscripts/PRJNA796958/SC.TotalCoverage50/ProteinsFiles/comp73852_c0_seq1.unique_proteins.csv.gz"
# firstcolpos = 16
# # fractions = [0.1]
# algs = ["Ascending", "Descending"]
# onlymaxdistinct = false
# maxmainthreads = 20
# fractions = [0.2, 0.4, 0.6, 0.8, 1.0]

# substitutionmatrix = nothing
# aagroups = nothing
# similarityscorecutoff = 0
# similarityvalidator = :(>=)


# distinctdf = prepare_distinctdf(
#     distinctfile, delim, innerdelim, truestrings, falsestrings
# )

# allprotsdf, firstcolpos = prepare_allprotsdf!(
#     allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos
# )

# # the possible amino acids each protein has in each position
# M = Matrix(allprotsdf[:, firstcolpos:end])
# # the distances between any two proteins according to `M`
# Δ = begin
#     if substitutionmatrix !== nothing
#         distances(M, substitutionmatrix, similarityscorecutoff, similarityvalidator)
#     elseif aagroups !== nothing
#         distances(M, aagroups)
#     else
#         distances(M)
#     end
# end

# # considering only desired solutions (rows' indices)
# solutions = choosesolutions(distinctdf, fractions, algs, onlymaxdistinct)

# # # allsubsolutions = collect(Iterators.partition(solutions, maxmainthreads))

# minmainthreads = minimum([Int(Threads.nthreads() / 5), Int(length(solutions) / 4)])
# allsubsolutions = collect(Iterators.partition(solutions, minmainthreads))

# # allprotsdf[:, firstcolpos-4:firstcolpos-1]
# # @assert all(length.(allprotsdf[!, "AdditionalSupportingReadsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])
# # @assert all(length.(allprotsdf[!, "AdditionalSupportingProteinsIDs"]) .== [0 for _ ∈ 1:size(allprotsdf, 1)])

# # # results_1_4 = additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, allsubsolutions[1])

# # # results = [
# # #     additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
# # #     for subsolutions ∈ allsubsolutions
# # # ]

# solution = 1


# solutionrow = distinctdf[solution, :]

# baseallprotsdf = deepcopy(allprotsdf[:, begin:firstcolpos-1])

# # if "AvailableReads" ∈ names(solutionrow) && solutionrow["Fraction"] < 1.0
# if "AvailableReads" ∈ names(solutionrow)
#     availablereads = solutionrow["AvailableReads"]
#     allreadsperprotein = baseallprotsdf[!, "Reads"]
#     availablereadsperprotein = [
#         [read for read ∈ reads if read ∈ availablereads]
#         for reads ∈ allreadsperprotein
#     ]
#     baseallprotsdf[:, "Reads"] .= availablereadsperprotein
#     baseallprotsdf[:, "NumOfReads"] .= length.(availablereadsperprotein)
#     baseallprotsdf = filter("NumOfReads" => x -> x > 0, baseallprotsdf)
# end





# one_result = one_solution_additional_assignment_considering_available_reads(distinctdf, allprotsdf, firstcolpos, Δ, solution)

# results = tcollect(
#     additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
#     for subsolutions ∈ allsubsolutions
# )
# finalresults = vcat((skipmissing(Iterators.flatten(results))...))






























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


