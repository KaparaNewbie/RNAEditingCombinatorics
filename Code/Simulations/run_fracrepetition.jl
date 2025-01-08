using DataFrames
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
    seed,
    algrepetitions::Int64,
    run_solve_threaded::Bool,
    sortresults::Bool,
    algs::Vector{String},
)
    @info "$(loggingtime())\trun_fracrepetition" fraction nsamplerows fracrepetition algrepetitions run_solve_threaded sortresults algs myid()
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
        sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol, randseed)
    else
        sampleG, availablereads = get_graph_sample_and_available_reads(G, fraction, nsamplerows, df, idcol)
    end

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