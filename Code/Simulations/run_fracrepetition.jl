using DataFrames
using CSV
using Transducers
import Base.Threads.@spawn
using Distributed


"Distributed version, w/o algrepetitions, so threads are only needed for subgraph sampling"
function run_fracrepetition(
	df::DataFrame,
	idcol::String,
	ArrG::DArray,
	fraction::Float64,
	nsamplerows::Int64,
	fracrepetition::Int64,
	consistentfracsampling::Bool,
	samplerows::Vector{Int64},
	# randseed::Int,
	algrepetitions::Int64,
	run_solve_threaded::Bool,
	sortresults::Bool,
	algs::Vector{String},
)
	@info "$(loggingtime())\trun_fracrepetition" fraction nsamplerows fracrepetition consistentfracsampling algrepetitions run_solve_threaded sortresults algs myid()
	# assemble sub neighborhood lists of indistinguishable sampled rows by using the pre-computed complete graph

	# print(df)

	# G = @timeit to "`G`" ArrG[1]  # retrive G which is the single element in the distributed array ArrG
	G = ArrG[1]  # retrive G which is the single element in the distributed array ArrG

	# G = try
	#     @timeit to "`G`" ArrG[1]
	# catch e
	#     @warn "$(loggingtime())\tcan't extract G from ArrG in run_fracrepetition" fraction fracrepetition ArrG e
	#     return
	# end

	# fraction = 0.6
	# nrows = size(df, 1)
	# nsamplerows = convert(Int, round(fraction * nrows))

	# sampleG = @timeit to "get_graph_sample" get_graph_sample(G, fraction, nsamplerows, df, idcol)
	# sampleG = get_graph_sample(G, fraction, nsamplerows, df, idcol)
	# while

	if consistentfracsampling
		# sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol, randseed)
		sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, samplerows, df, idcol)
	else
		sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol)
	end
	# @info typeof(availablereads)

	# sampleG = try
	#     @timeit to "get_graph_sample" get_graph_sample(G, fraction, nsamplerows, df, idcol)
	# catch e
	#     @warn "$(loggingtime())\tcan't extract sampleG from G in run_fracrepetition" fraction fracrepetition e
	#     return
	# end

	# obtain sets of distinct rows
	# results = @timeit to "solve" solve(
	results = solve(
		sampleG,
		fraction,
		fracrepetition,
		algrepetitions,
		run_solve_threaded,
		sortresults,
		algs,
	)

	results[!, "AvailableReads"] .= join(availablereads, ",")

	# results = try
	#     @timeit to "solve" solve(
	#         sampleG,
	#         fraction,
	#         fracrepetition,
	#         algrepetitions,
	#         run_solve_threaded,
	#         sortresults,
	#         algs,
	#     )
	# catch e
	#     @warn "$(loggingtime())\tcan't solve sampleG in run_fracrepetition" fraction fracrepetition e
	#     return
	# end


	# garbage collection
	sampleG = nothing # todo check if not assigning nothing to sampleG helps prevent "ArgumentError: array must be non-empty"
	GC.gc()
	return results
end







"A test version of run_fracrepetition which writes the names of sampled reads to a csv file - used for debugging of consistentfracsampling."
function run_fracrepetition(
	df::DataFrame,
	idcol::String,
	ArrG::DArray,
	fraction::Float64,
	nsamplerows::Int64,
	fracrepetition::Int64,
	consistentfracsampling::Bool,
	samplerows::Vector{Int64},
	algrepetitions::Int64,
	run_solve_threaded::Bool,
	sortresults::Bool,
	algs::Vector{String},
	outdir::String,
	samplename::String,
)
	# @info "$(loggingtime())\trun_fracrepetition" fraction nsamplerows fracrepetition consistentfracsampling randseed algrepetitions run_solve_threaded sortresults algs outdir postfix_to_add myid()
	@info "$(loggingtime())\trun_fracrepetition" fraction fracrepetition nsamplerows consistentfracsampling algrepetitions run_solve_threaded sortresults algs outdir samplename myid()
	# assemble sub neighborhood lists of indistinguishable sampled rows by using the pre-computed complete graph

	# print(df)

	# G = @timeit to "`G`" ArrG[1]  # retrive G which is the single element in the distributed array ArrG
	G = ArrG[1]  # retrive G which is the single element in the distributed array ArrG

	if consistentfracsampling
		# sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol, randseed)
		sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, samplerows, df, idcol)
	else
		sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol)
	end

	# write the names of sampled reads to a csv file

	# print(typeof(availablereads))
	@info typeof(availablereads)

	# sampledreadsoutfile = joinpath(abspath(outdir), "$samplename.SampledReads.$fraction.$postfix_to_add.csv")
	sampledreadsoutfile = joinpath(abspath(outdir), "$samplename.SampledReads.$fraction.csv")
	CSV.write(sampledreadsoutfile, DataFrame(:Read => availablereads))

	sampledrowsoutfile = joinpath(abspath(outdir), "$samplename.SampledRows.$fraction.csv")
	CSV.write(sampledrowsoutfile, DataFrame(:Row => samplerows))

	# obtain sets of distinct rows
	# results = @timeit to "solve" solve(
	results = solve(
		sampleG,
		fraction,
		fracrepetition,
		algrepetitions,
		run_solve_threaded,
		sortresults,
		algs,
	)

	results[!, "AvailableReads"] .= join(availablereads, ",")



	# garbage collection
	sampleG = nothing # todo check if not assigning nothing to sampleG helps prevent "ArgumentError: array must be non-empty"
	GC.gc()
	return results
end