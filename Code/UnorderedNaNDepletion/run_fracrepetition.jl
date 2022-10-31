using DataFrames
using Transducers
import Base.Threads.@spawn
using Distributed


"Distributed version, w/o algrepetitions, so threads are only needed for subgraph sampling"
function run_fracrepetition(
    df::DataFrame,
    idcol::String,
    # G::Dict,
    ArrG::DArray,
    fraction::Float64,
    nsamplerows::Int64,
    fracrepetition::Int64,
    algrepetitions::Int64,
    run_solve_threaded::Bool,
    sortresults::Bool,
    algs::Vector{String},
    # algfuncs::Dict{String,Function}
)
    @info "$(loggingtime())\trun_fracrepetition" fraction nsamplerows fracrepetition algrepetitions run_solve_threaded sortresults algs myid()
    # assemble sub neighborhood lists of indistinguishable sampled rows by using the pre-computed complete graph
    G = @timeit to "`G`" ArrG[1]  # retrive G which is the single element in the distributed array ArrG
    sampleG = @timeit to "get_graph_sample" get_graph_sample(G, fraction, nsamplerows, df, idcol)
    # obtain sets of distinct rows
    # results = solve(sampleG, fraction, fracrepetition)
    results = @timeit to "solve" solve(
        sampleG,
        fraction,
        fracrepetition,
        algrepetitions,
        run_solve_threaded,
        sortresults,
        algs,
        # algfuncs
    )
    # garbage collection
    sampleG = nothing
    GC.gc()
    return results
end