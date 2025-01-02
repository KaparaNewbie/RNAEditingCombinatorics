using ArgParse
using BioSequences # for BioSequences.toAAset
using BioAlignments  # right now it's a local branch from my own repo - for SubstitutionMatrix
using DataFrames
using CSV
using Random # for Rand.seed
using StatsBase  # for StatsBase.sample
using InlineStrings  # for saving space on `idcol` (Read/Protein) in input `df`
using ThreadsX
using Transducers


include(joinpath(@__DIR__, "consts.jl")) # for ∅ & AA_groups
include(joinpath(@__DIR__, "timeformatters.jl"))
include(joinpath(@__DIR__, "preparedf.jl"))
include(joinpath(@__DIR__, "indistinguishable_rows.jl"))

"""
Define command-line arguments.
"""
function parsecmd()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--complete_infiles"
		help = "One or more csv files representing unique proteins of a ground-truth dataset."
		nargs = '+'
		action = :store_arg
		required = true
		"--errored_na_infiles"
		help = "Corresponding csv files representing unique proteins of coupled test datasets."
		nargs = '+'
		action = :store_arg
		required = true
		"--out_files"
		help = "Corresponding output files. Will undergo gzip compression if the file extension is `.gz`."
		nargs = '+'
		action = :store_arg
		required = true
		"--delim"
		help = "Delimiter for input/output csv files."
		arg_type = String
		default = "\t"
		# "--prefix_to_remove", "--prefix"
		# help = "Remove `prefix` from output files` names."
		# default = ""
		# "--postfix_to_remove", "--postfix"
		# help = "Remove `postfix` from output files` names, e.g., `.unique_reads.csv`."
		# default = ""
		# "--postfix_to_add"
		# help = "Add `postfix` to output files' names, e.g., `\$sample.DistinctUnique{Reads,Proteins}\$postfix.\$time.csv`."
		# default = ""
		"--idcol"
		help = "Label of unique samples columns. Typically `Transcripts` for `Reads` and `Protein` for `Proteins`."
		default = "Protein"
		# required = true
		"--firstcolpos"
		help = "Int location of the first editing position column of each file in `infiles`. As of now, should be 9 for `Reads` and 15 for `Proteins`."
		arg_type = Int
		default = 15
		# required = true
		"--testfraction"
		help = "Fraction of the dataset to be used for testing. That fraction will correspond to `maxfrac == 1.0`."
		arg_type = Float64
		default = 1.0
		range_tester = x -> 0.0 < x <= 1.0
		"--randseed"
		help = "Random seed for sampling test data."
		arg_type = Int
		default = 1892
	end
	return parse_args(s)
end


# complete_infile = "/private7/projects/Combinatorics/Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.Complete.UniqueProteins.tsv.gz"
# errored_na_infile = "/private7/projects/Combinatorics/Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.Errored.PartiallyUnknown.UniqueProteins.tsv.gz"
# out_file = "/private7/projects/Combinatorics/Simulations/GraphAssessment/comp141434_c0_seq1.CBPC1_HUMAN.UP0_09.Rep1.FalsePositives.tsv.gz"

# delim = "\t"
# # datatype = "Proteins"
# idcol = "Protein"
# firstcolpos = 15
# testfraction = 0.1 # for testing purposes only!
# # testfraction = 1.0
# # substitutionmatrix = nothing
# # aagroups = nothing


# # randseed = 1892
# complete_protein_suffix = "-C"
# errored_na_protein_suffix = "-D"


function preparecoupleddf!(
	complete_infile::String, errored_na_infile::String, delim::String, firstcolpos::Int,
	testfraction::Float64 = 1.0, idcol::String = "Protein", complete_protein_suffix::String = "-C", errored_na_protein_suffix::String = "-D", randseed = 1892,
)
	@info "$(loggingtime())\tpreparecoupleddf!" complete_infile errored_na_infile delim firstcolpos testfraction idcol complete_protein_suffix errored_na_protein_suffix randseed

	# read the files into a df

	infiles = [complete_infile, errored_na_infile]
	protein_suffixes = [complete_protein_suffix, errored_na_protein_suffix]

	dfs = []

	# infile = infiles[1]
	# protein_suffix = protein_suffixes[1]

	for (infile, protein_suffix) in zip(infiles, protein_suffixes)
		# note: "Protein" is actually the ID of a *unique* protein!!!
		df1 = DataFrame(CSV.File(infile, delim = delim, select = collect(1:firstcolpos-1), types = Dict("Protein" => String, "Reads" => String)))
		# df1[!, "Protein"] = InlineString.(df1[!, :Protein])
		df1[!, "Protein"] = InlineString.(df1[!, :Protein] .* protein_suffix)
		# make sure columns of AAs containing only Ts aren't parsed as boolean columns
		df2 = DataFrame(CSV.File(infile, delim = delim, drop = collect(1:firstcolpos-1), types = String))
		df = hcat(df1, df2)


		# take a subset of the df for testing purposes
		if testfraction < 1.0
			# Random.seed!(randseed)
			nrows = size(df, 1)
			nsamplerows = Int(round(testfraction * nrows))
			# samplerows = sample(1:nrows, nsamplerows, replace=false)
			samplerows = sample(MersenneTwister(randseed), 1:nrows, nsamplerows, replace = false)
			df = df[samplerows, :]
		end

		# flatten the df by exploding the `Reads` col, 
		# which denotes the reads supporting the unique observation (read / unique) the row represents
		transform!(df, :Reads => (x -> split.(x, ",")) => :Reads)
		df = flatten(df, :Reads)
		# include also the name of every read supporting a unique protein
		rename!(df, "Reads" => "Read")
		df = hcat(select(df, idcol, "Read"), toAAset.(df[:, firstcolpos:end]))

		push!(dfs, df)
	end

	# dfs[1]
	# dfs[2]
	@assert length(dfs) == 2

	firstcolpos = 3

	coupleddf = vcat(dfs..., cols = :union)


	# remove uniformative cols 
	# (altough we also apply the function on the idcol it shouldn't matter for the idcol, 
	# if it has more than one unique value)
	# informativedf = df[:, map(col -> length(unique(col)) > 1, eachcol(df))]
	# if size(informativedf)[1] == 0
	#     df = df[1, :]
	# end
	coupleddf = coupleddf[:, map(col -> length(unique(col)) > 1, eachcol(coupleddf))]

	return coupleddf, firstcolpos
end


function coupled_indistinguishable_rows(coupleddf, idcol, firstcolpos)

	@info "$(loggingtime())\tcoupled_indistinguishable_rows"

	M, ids = preprocess(coupleddf, idcol, firstcolpos)

	nrows = size(M, 1)
	length(ids) == length(ThreadsX.unique(ids)) == nrows || error(
		"Each id (a vertex in the graph) should be unique. " *
		"There should be a bijection between M's rows and ids, in order of appearence.",
	)
	nodeseltype = eltype(ids)

	AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
	# sets whose intersection is empty
	distinctAAsets = ThreadsX.Set(
		[
		(x, y)
		for x ∈ AAsets
		for y ∈ AAsets
		if x ∩ y == ∅
	]
	)

	# Gs = tcollect(indistinguishable_vs_for_u(M, distinctAAsets, ids, nodeseltype, i) for i ∈ 1:nrows)
	Gs = Vector{Any}(undef, nrows)
	Threads.@threads for i ∈ 1:nrows
		Gs[i] = indistinguishable_vs_for_u(M, distinctAAsets, ids, nodeseltype, i)
	end
	G = Dict(k => v for g ∈ Gs for (k, v) ∈ g) # merging Gs into G; each g has a single (k, v) pair
	return G
end


function removesuffix(str, suffix)
	if !occursin(suffix, str)
		return str
	end
	return str[1:end-length(suffix)]
end


function custom_reduce(strings)
	if length(strings) == 0
		return ""
	else
		return reduce((x, y) -> x * "," * y, strings)
	end
end


function run_sample(
	complete_infile, errored_na_infile, out_file, delim, firstcolpos,
	testfraction, idcol, complete_protein_suffix, errored_na_protein_suffix, randseed,
)
	@info "$(loggingtime())\trun_sample" complete_infile errored_na_infile out_file delim firstcolpos testfraction idcol complete_protein_suffix errored_na_protein_suffix randseed

	coupleddf, firstcolpos = preparecoupleddf!(
		complete_infile, errored_na_infile, delim, firstcolpos,
		testfraction, idcol, complete_protein_suffix, errored_na_protein_suffix, randseed,
	)

	G_coupled = try
		coupled_indistinguishable_rows(coupleddf, idcol, firstcolpos)
	catch e
		@warn "$(loggingtime())\tcoupled_indistinguishable_rows failed for $complete_infile and $errored_na_infile" e
		return
	end

	G_D = Dict(k => v for (k, v) in G_coupled if occursin(errored_na_protein_suffix, k))

	G_FP = Dict()
	# for each damaged protein
	for k in keys(G_D)
		vs = []
		# add its indistinguishable complete proteins to the FP dict
		for v in values(G_D[k])
			if occursin(complete_protein_suffix, v)
				# remove the suffix from the complete protein
				v = removesuffix(v, complete_protein_suffix)
				push!(vs, v)
			end
		end
		# remove the suffix from the damaged protein
		k = removesuffix(k, errored_na_protein_suffix)
		G_FP[k] = vs
	end

	ks = [k for k in keys(G_FP)]
	n_neighbours = length.(values(G_FP))
	vs = [custom_reduce(v) for v in values(G_FP)]

	df = DataFrame(
		"Errored+PartiallyUnknownProtein" => ks,
		"NumOfIndistinguishableCompleteProteins" => n_neighbours,
		"IndistinguishableCompleteProteins" => vs,
	)
	sort!(df, "Errored+PartiallyUnknownProtein")

	compress = occursin(".gz", out_file) ? true : false

	CSV.write(out_file, df; delim, compress)

end



function main()

	# read command-line args
	parsedargs = parsecmd()

	complete_infiles = parsedargs["complete_infiles"]
	errored_na_infiles = parsedargs["errored_na_infiles"]
	out_files = parsedargs["out_files"]
	delim = parsedargs["delim"]
	# prefix_to_remove = parsedargs["prefix_to_remove"]
	# postfix_to_remove = parsedargs["postfix_to_remove"]
	# postfix_to_add = parsedargs["postfix_to_add"]
	idcol = parsedargs["idcol"]
	firstcolpos = parsedargs["firstcolpos"]
	testfraction = parsedargs["testfraction"]
	randseed = parsedargs["randseed"]

	# @info "$(loggingtime())\tmain" complete_infiles errored_na_infiles out_files delim idcol firstcolpos testfraction randseed
	@assert length(complete_infiles) == length(errored_na_infiles) == length(out_files)

	num_of_datasets = length(complete_infiles)

	@info "$(loggingtime())\tmain" num_of_datasets idcol firstcolpos testfraction randseed

	complete_protein_suffix = "-C"
	errored_na_protein_suffix = "-D"

	for (complete_infile, errored_na_infile, out_file) in zip(complete_infiles, errored_na_infiles, out_files)
		run_sample(
			complete_infile, errored_na_infile, out_file, delim, firstcolpos,
			testfraction, idcol, complete_protein_suffix, errored_na_protein_suffix, randseed,
		)
	end

end


main()


# coupleddf, firstcolpos = preparecoupleddf!(
# 	complete_infile, errored_na_infile, delim, firstcolpos,
# 	testfraction,
# )

# # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
# G_coupled = try
# 	coupled_indistinguishable_rows(coupleddf, idcol, firstcolpos)
# catch e
# 	@warn "$(loggingtime())\tcoupled_indistinguishable_rows failed for $complete_infile and $errored_na_infile" e
# 	return
# end


# G_D = Dict(k => v for (k, v) in G_coupled if occursin(errored_na_protein_suffix, k))


# G_FP = Dict()
# # for each damaged protein
# for k in keys(G_D)
# 	vs = []
# 	# add its indistinguishable complete proteins to the FP dict
# 	for v in values(G_D[k])
# 		if occursin(complete_protein_suffix, v)
# 			# remove the suffix from the complete protein
# 			v = removesuffix(v, complete_protein_suffix)
# 			push!(vs, v)
# 		end
# 	end
# 	# remove the suffix from the damaged protein
# 	k = removesuffix(k, errored_na_protein_suffix)
# 	G_FP[k] = vs
# end



# ks = [k for k in keys(G_FP)]
# n_neighbours = length.(values(G_FP))
# vs = [custom_reduce(v) for v in values(G_FP)]

# df = DataFrame(
# 	"Errored+PartiallyUnknownProtein" => ks,
# 	"NumOfIndistinguishableCompleteProteins" => n_neighbours,
# 	"IndistinguishableCompleteProteins" => vs,
# )
# sort!(df, "Errored+PartiallyUnknownProtein")

# compress = occursin(".gz", out_file) ? true : false
# CSV.write(out_file, df; delim, compress)








# using DataStructures
# using CairoMakie
# CairoMakie.activate!()


# # fig = Figure(; size = (1200, 800), fonts = (; regular= "sans"), fontsize = 20)

# fig = Figure()
# ax1 = Axis(
# 	fig[1, 1];
# 	xlabel = "# of indistinguishable complete proteins\nper errored + partially-unknown protein",
# 	ylabel = "# of errored + partially-unknown proteins",
# )
# hist!(
# 	ax1,
# 	n_neighbours,
# 	bins = 20, color = :red, strokewidth = 1, strokecolor = :black,
# )
# fig



# fp_counter = counter(length.(values(G_FP)))
# xs = []
# ys = []
# for (k, v) in fp_counter
# 	push!(xs, k)
# 	push!(ys, v)
# end

# fig = Figure()
# ax1 = Axis(
# 	fig[1, 1];
# 	xlabel = "# of indistinguishable complete proteins\nper errored + partially-unknown protein",
# 	ylabel = "# of errored + partially-unknown proteins",
# )
# barplot!(
# 	ax1,
# 	xs, ys,
# 	color = :red, strokewidth = 1, strokecolor = :black,
# )
# fig
