# todo cloud change, verify on the server as well
n_workers = length(workers())
if n_workers > 1
    # remove previous workers
    rmprocs(workers())
    # add new one, but with limited number of threads
    addprocs(
        n_workers,
        exeflags=[
            "--threads=$(Int(round(Threads.nthreads() / n_workers)))",
            "--project"
        ]
    )
end


# https://discourse.julialang.org/t/julia-python-equivalent-of-main/35433
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    # Pkg.activate(".") # todo cloud change, verify on the server as well
else
    using BenchmarkTools
end


using ArgParse
using DataFrames # for reading input and writing output
using CSV
using DelimitedFiles
using LinearAlgebra
using StatsBase  # for StatsBase.sample
using IterTools  # for IterTools.imap
using Random # for MersenneTwister & shuffle
using BioSequences # for BioSequences.toAAset
using Dates # for logging
import Base.Threads.@spawn
using Distributed
using Transducers
using ThreadsX
using InlineStrings  # for saving space on "Reads" col in input `df`


# todo cloud change, verify on the server as well
# n_workers = length(workers())
# if n_workers > 1
#     # remove previous workers
#     rmprocs(workers())
#     # add new one, but with limited number of threads
#     addprocs(
#         n_workers,
#         exeflags=[
#             "--threads=$(Int(round(Threads.nthreads() / n_workers)))",
#             "--project"
#         ]
#     )
# end

include(joinpath(@__DIR__, "consts.jl")) # for ∅
include(joinpath(@__DIR__, "timeformatters.jl"))
include(joinpath(@__DIR__, "preparedf.jl"))
include(joinpath(@__DIR__, "indistinguishable_rows.jl"))


@everywhere begin
    using DataFrames # for reading input and writing output
    using Distributed, DistributedArrays
    using Transducers
    using ThreadsX
    using InlineStrings  # got an error without it
    using Dates # for logging
    using StatsBase  # for StatsBase.sample
    using IterTools  # for IterTools.imap
    using Random # for MersenneTwister & shuffle
    include(joinpath(@__DIR__, "consts.jl")) # for ∅
    include(joinpath(@__DIR__, "timeformatters.jl"))
    include(joinpath(@__DIR__, "solve.jl"))
    include(joinpath(@__DIR__, "sample.jl"))
    include(joinpath(@__DIR__, "run_fracrepetition.jl"))
    # using Logging, LoggingExtras
    # include(joinpath(@__DIR__, "setlogger.jl"))
end


"""
Define command-line arguments.
"""
function parsecmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--infiles"
        help = "One or more csv files representing unique reads/proteins."
        nargs = '+'
        action = :store_arg
        required = true
        "--delim"
        help = "Delimiter for input/output csv files."
        arg_type = String
        default = "\t"
        "--prefix_to_remove", "--prefix"
        help = "Remove `prefix` from output files` names."
        default = ""
        "--postfix_to_remove", "--postfix"
        help = "Remove `postfix` from output files` names, e.g., `.unique_reads.csv`."
        default = ""
        "--postfix_to_add"
        help = "Add `postfix` to output files` names, e.g., `\$sample.DistinctUnique{Reads,Proteins}\$postfix.\$time.csv`."
        default = ""
        "--idcol"
        help = "Label of unique samples columns. Typically `Transcripts` for `Reads` and `Protein` for `Proteins`."
        required = true
        "--firstcolpos"
        help = "Int location of the first editing position column of each file in `infiles`. As of now, should be 9 for `Reads` and 15 for `Proteins`."
        arg_type = Int
        required = true
        "--datatype"
        help = "Data type of the input files. Either `Reads` or `Proteins`."
        required = true
        "--outdir"
        help = "Write output files to this directory."
        required = true
        "--fracstep"
        help = "Step size for the fraction of the dataset to be used for each iteration, s.t. `fractions = collect(fracstep:fracstep:1.0) ∪ 1.0`."
        arg_type = Float64
        default = 0.1
        "--maxfrac"
        help = "Maximum fraction of the dataset to be sampled."
        arg_type = Float64
        default = 1.0
        "--fracrepetitions"
        help = "Number of repetitions of the sampling procedure for each fraction of the data."
        arg_type = Int
        default = 10
        "--algrepetitions"
        help = "Number of repetitions for each algorithm is run on each repetition of each fraction of the data."
        arg_type = Int
        default = 5
        "--testfraction"
        help = "Fraction of the dataset to be used for testing. That fraction will correspond to `maxfrac == 1.0`."
        arg_type = Float64
        default = 1.0
        range_tester = x -> 0.0 < x <= 1.0
        "--randseed"
        help = "Random seed for sampling test data."
        arg_type = Int
        default = 1892
        "--run_solve_threaded"
        help = "Run different algorithms/algrepetitions in parallel using threads"
        action = :store_true
        "--sortresults"
        help = "Sort distinct unique samples of reads/proteins."
        action = :store_true
        "--algs"
        help = "Use these algorithms to obtain distinct unique samples of reads/proteins."
        default = ["Ascending", "Descending", "Unordered"]
        nargs = '+'
        range_tester = x -> x ∈ ["Ascending", "Descending", "Unordered"]
        "--gcp"
        help = "Program is run on a google cloud VM."
        action = :store_true
        "--shutdowngcp"
        help = "Shutdown google cloud VM when the program ends."
        action = :store_true

    end
    return parse_args(s)
end




"""
The final result: 10 repetions of `f1` and `f2` over the neighborhood matrix of the graph, for each fraction of reads.  
Each such result contains a list of compatible unique reads (they are all different from each other).
"""
function main()

    # read command-line args
    parsedargs = parsecmd()
    infiles = parsedargs["infiles"]
    delim = parsedargs["delim"]
    prefix_to_remove = parsedargs["prefix_to_remove"]
    postfix_to_remove = parsedargs["postfix_to_remove"]
    postfix_to_add = parsedargs["postfix_to_add"]
    idcol = parsedargs["idcol"]
    firstcolpos = parsedargs["firstcolpos"]
    datatype = parsedargs["datatype"]
    outdir = parsedargs["outdir"]
    fracstep = parsedargs["fracstep"]
    maxfrac = parsedargs["maxfrac"]
    fracrepetitions = parsedargs["fracrepetitions"]
    algrepetitions = parsedargs["algrepetitions"]
    testfraction = parsedargs["testfraction"]
    randseed = parsedargs["randseed"]
    run_solve_threaded = parsedargs["run_solve_threaded"]
    sortresults = parsedargs["sortresults"]
    algs = parsedargs["algs"]
    gcp = parsedargs["gcp"]
    shutdowngcp = parsedargs["shutdowngcp"]

    algs::Vector{String} = String.(algs)

    # infiles = ["D.pealeii/MpileupAndTranscripts/IlluminaExample/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv"]
    # prefix = "reads.sorted.aligned.filtered."
    # postfix = ".unique_proteins.csv"
    # delim = "\t"
    # idcol = "Protein"
    # editingcol = ""
    # firstcolpos = 15
    # removeuneditedrows = false
    # allrowsedited = false
    # datatype = "Proteins"
    # outdir = "D.pealeii/MpileupAndTranscripts/IlluminaExample"
    # fracstep = 0.5
    # maxfrac = 1.0
    # fracrepetitions = 5
    # algrepetitions = 2
    # testfraction = 0.01
    # run_fraction_parallelism = "sequential"
    # run_fracrepetition_parallelism = "distributed"
    # randseed = 1892

    @info "$(loggingtime())\tmain" delim prefix_to_remove postfix_to_remove postfix_to_add idcol firstcolpos datatype outdir fracstep maxfrac fracrepetitions algrepetitions testfraction randseed run_solve_threaded sortresults algs

    # do the thing for each file, in ascending order of file size
    sort!(infiles, by=(x -> filesize(x)))
    samplesnames = map(x -> replace(splitdir(x)[2], prefix_to_remove => "", postfix_to_remove => ""), infiles)

    # logfile = joinpath(outdir, "log.$(writingtime()).txt")
    # run(`touch $logfile`)
    # @everywhere setlogger($logfile) # https://discourse.julialang.org/t/copying-variable-to-all-remote-processes/26518/4?u=kaparanewbie

    for x ∈ eachindex(infiles)  # could be eachindex(samplesnames) as well
        infile = infiles[x]
        samplename = samplesnames[x]
        run_sample(
            infile, delim, samplename, idcol, firstcolpos, datatype, outdir, postfix_to_add,
            fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed,
            run_solve_threaded, sortresults, algs
        )
    end

    # shutdown gcp vm
    gcp && shutdowngcp && run(`sudo shutdown`)
end



"""
# Examples
```julia-repl
julia> interleaved_fold([1, 2, 3, 4])
[1, 4, 2, 3]
```
```julia-repl
julia> interleaved_fold([1, 2, 3, 4, 5])
[1, 5, 2, 4, 3]
```
"""
function interleaved_fold(A)
    middle = Int(ceil(length(A) / 2))
    if iseven(length(A))
        a = A[begin:middle]
        b = reverse(A[middle+1:end])
        Y = collect(Iterators.flatten(zip(a, b)))
    else # odd
        a = A[begin:middle-1]
        b = reverse(A[middle+1:end])
        Y = vcat(collect(Iterators.flatten(zip(a, b))), A[middle])
    end
    Y
end


"""
For `fracrepetitions = 2`, `fracstep = 0.3`, `maxfrac = 1.0`, 
and `df` with `length(df) == 100`, the output is:

(fraction = 0.25, nsamplerows = 25, fracrepetition = 1)  
 (fraction = 1.0, nsamplerows = 100, fracrepetition = 1)  
 (fraction = 0.5, nsamplerows = 50, fracrepetition = 1)  
 (fraction = 0.75, nsamplerows = 75, fracrepetition = 1)  
 (fraction = 0.25, nsamplerows = 25, fracrepetition = 2)  
 (fraction = 1.0, nsamplerows = 100, fracrepetition = 2)  
 (fraction = 0.5, nsamplerows = 50, fracrepetition = 2)  
 (fraction = 0.75, nsamplerows = 75, fracrepetition = 2)  
"""
function prep_pmap_input(fracstep, maxfrac, df, fracrepetitions)
    fractions = collect(fracstep:fracstep:maxfrac) ∪ maxfrac
    # define nsamplerows for each fraction
    nrows = size(df, 1)
    fraction_nsamplerows = [convert(Int, round(fraction * nrows)) for fraction ∈ fractions]

    fracrepetitions_inputs = [
        interleaved_fold(
            [
            (fraction=fraction, nsamplerows=nsamplerows, fracrepetition=fracrepetition)
            for (fraction, nsamplerows) ∈ zip(fractions, fraction_nsamplerows)
        ]
        )
        for fracrepetition ∈ 1:fracrepetitions
    ]
    fracrepetitions_inputs = vcat(fracrepetitions_inputs...)

    return fracrepetitions_inputs
end


# fracstep = 0.25
# maxfrac = 1.0
# nrows = 100
# fracrepetitions = 2
# fractions = collect(fracstep:fracstep:maxfrac) ∪ maxfrac
# fraction_nsamplerows = [convert(Int, round(fraction * nrows)) for fraction ∈ fractions]
# fracrepetitions_inputs = [
#     interleaved_fold(
#         [
#         (fraction=fraction, nsamplerows=nsamplerows, fracrepetition=fracrepetition)
#         for (fraction, nsamplerows) ∈ zip(fractions, fraction_nsamplerows)
#     ]
#     )
#     for fracrepetition ∈ 1:fracrepetitions
# ]
# fracrepetitions_inputs = vcat(fracrepetitions_inputs...)




function run_sample(
    infile::String,
    delim::String,
    samplename::String,
    idcol::String,
    firstcolpos::Int,
    datatype::String,
    outdir::String,
    postfix_to_add::String,
    fracstep::Float64,
    maxfrac::Float64,
    fracrepetitions::Int,
    algrepetitions::Int,
    testfraction::Float64,
    randseed::Int,
    run_solve_threaded::Bool,
    sortresults::Bool,
    algs::Vector{String}
)
    @info "$(loggingtime())\trun_sample" infile samplename

    df, firstcolpos = try
        preparedf!(
            infile, delim, datatype, idcol, firstcolpos,
            testfraction, randseed
        )
    catch e
        @warn "$(loggingtime())\tpreparedf! failed for $infile" e
        return
    end

    # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
    G = try
        indistinguishable_rows(df, idcol; firstcolpos)
    catch e
        @warn "$(loggingtime())\tindistinguishable_rows failed for $infile" e
        return
    end
    ArrG = @DArray [G for _ ∈ 1:1]  # move G across processes on a distributed array in order to save memory

    # having built G, we only need to keep the reads and unique reads/proteins they support
    select!(df, idcol)
    GC.gc() # free up memory, just in case


    # define input parameters for pmap
    fracrepetitions_inputs = prep_pmap_input(fracstep, maxfrac, df, fracrepetitions)

    fracrepetitions_results = pmap(
        fracrepetitions_inputs;
        retry_delays=ExponentialBackOff(n=3, first_delay=5, max_delay=1000)
    ) do (fraction, nsamplerows, fracrepetition)
        try
            # run_fracrepetition(df, idcol, G, fraction, fracrepetition, nsamplerows)
            run_fracrepetition(
                df,
                idcol,
                ArrG,
                fraction,
                nsamplerows,
                fracrepetition,
                algrepetitions,
                run_solve_threaded,
                sortresults,
                algs
            )

            # run_fracrepetition(
            #     df::DataFrame,
            #     idcol::String,
            #     # G::Dict,
            #     ArrG::DArray,
            #     fraction::Float64,
            #     nsamplerows::Int64,
            #     fracrepetition::Int64,
            #     algrepetitions::Int64,
            #     run_solve_threaded::Bool,
            #     sortresults::Bool,
            #     algs::Vector{String},
            #     # algfuncs::Dict{String,Function}
            # )

        catch e
            @warn "$(loggingtime())\trun_fracrepetition failed" fraction fracrepetition e
            missing
        end
    end

    results::DataFrame = vcat(skipmissing(fracrepetitions_results)...)
    sort!(results, ["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])

    # write the results to a csv file
    # outfile = abspath(outdir) * "/" * samplename * ".DistinctUnique$datatype.$(writingtime()).csv"
    outfile = joinpath(abspath(outdir), "$samplename.DistinctUnique$datatype$postfix_to_add.$(writingtime()).csv")
    CSV.write(outfile, results; delim)
end



main()



# @benchmark main()


# idcol = "Protein"
# delim = "\t"
# datatype = "Proteins"
# firstcolpos = 15
# fraction = 1.0
# testfraction = 1.0


# newG = let
#     infile = "D.pealeii/MpileupAndTranscripts/IlluminaExample2/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv"
#     firstcolpos = 15
#     outdir = "D.pealeii/MpileupAndTranscripts/IlluminaExample2"
#     df, firstcolpos = preparedf!(infile, delim, datatype, idcol, firstcolpos, testfraction)
#     # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
#     G = indistinguishable_rows(df, idcol; firstcolpos)
# end

# newdegsdf = let G = newG
#     vertices = Vector{eltype(keys(G))}(undef, length(G))
#     degrees = Vector{Int}(undef, length(G))
#     for (x, (k, v)) ∈ enumerate(G)
#         vertices[x] = k
#         degrees[x] = length(v)
#     end
#     degsdf = DataFrame("Node" => vertices, "Neighbors" => degrees)
#     sort!(degsdf, "Neighbors")
# end

# varinfo(r"newG") # 14.443 GiB Dict{String3, Set{String3}} with 170080 entries

# # CSV.write("/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/IlluminaExample2/newdegsdf.csv", newdegsdf)

# oldG = let
#     infile = "D.pealeii/MpileupAndTranscripts/IlluminaExample/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv"
#     firstcolpos = 15
#     outdir = "D.pealeii/MpileupAndTranscripts/IlluminaExample"
#     df, firstcolpos = preparedf!(infile, delim, datatype, idcol, firstcolpos, testfraction)
#     # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
#     G = indistinguishable_rows(df, idcol; firstcolpos)
# end


# olddegsdf = let G = oldG
#     vertices = Vector{eltype(keys(G))}(undef, length(G))
#     degrees = Vector{Int}(undef, length(G))
#     for (x, (k, v)) ∈ enumerate(G)
#         vertices[x] = k
#         degrees[x] = length(v)
#     end
#     degsdf = DataFrame("Node" => vertices, "Neighbors" => degrees)
#     sort!(degsdf, "Neighbors")
# end




# # CSV.write("/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/IlluminaExample/olddegsdf.csv", olddegsdf)

# varinfo(r"oldG") # oldG 110.567 GiB Dict{String31, Set{String31}} with 170966 entries


# n1 = olddegsdf[!, :Neighbors]
# n2 = newdegsdf[!, :Neighbors]

# minimum(n1)
# minimum(n2)
# mean(n1)
# mean(n2)
# median(n1)
# median(n2)
# maximum(n1)
# maximum(n2)

# s1 = sum(n1)
# s2 = sum(n2)
# min(s1, s2) / max(s1, s2)

# vertices = Vector{eltype(keys(G))}(undef, length(G))
# degrees = Vector{Int}(undef, length(G))
# for (x, (k, v)) ∈ enumerate(G)
#     vertices[x] = k
#     degrees[x] = length(v)
# end
# # degs_df = DataFrame("Node" => vertices, "Neighbors" => degrees)

# G_size = 16.576 # GiB varinfo(r"G")

# vs_plus_ns = length(vertices) + sum(degrees)
# bytes_needed = vs_plus_ns * 4
# gb_needed = bytes_needed / 10^9



# infile = "D.pealeii/MpileupAndTranscripts/IlluminaExample2/reads.sorted.aligned.filtered.comp141044_c0_seq2.unique_proteins.csv"
# samplename = "comp141044_c0_seq2"
# idcol = "Protein"
# outdir = "D.pealeii/MpileupAndTranscripts/IlluminaExample2"
# fracstep = 1.0
# maxfrac = 1.0
# fracrepetitions = 1
# algrepetitions = 1
# delim = "\t"
# datatype = "Proteins"
# fraction = 1.0
# fracrepetition = 1
# algrepetitions = 10
# firstcolpos = 15
# testfraction = 1.0

# df, firstcolpos = preparedf!(infile, delim, datatype, idcol, firstcolpos, testfraction)

# # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
# @time G = indistinguishable_rows(df, idcol; firstcolpos)

# vertices = Vector{eltype(keys(G))}(undef, length(G))
# degrees = Vector{Int}(undef, length(G))
# for (x, (k, v)) ∈ enumerate(G)
#     vertices[x] = k
#     degrees[x] = length(v)
# end
# # degs_df = DataFrame("Node" => vertices, "Neighbors" => degrees)

# G_size = 16.576 # GiB varinfo(r"G")

# vs_plus_ns = length(vertices) + sum(degrees)
# bytes_needed = vs_plus_ns * 4
# gb_needed = bytes_needed / 10^9

# # # having built G, we only need to keep the reads and unique reads/proteins they support
# # select!(df, idcol)

# # free()
# # GC.gc()
# # free()

# # nrows = size(df, 1)
# # nsamplerows = convert(Int, round(fraction * nrows)) # numbe

# # # assemble sub neighborhood lists of uncompatible sampled rows by using the pre-computed complete graph
# # sampleG = get_graph_sample(G, fraction, nsamplerows, df, idcol)

# # # # obtain $algrepetitions × 2 (the asc- and desc algorithms) results of compatible rows
# # @time results = solve(sampleG, fraction, fracrepetition, algrepetitions)
# # @benchmark solve(sampleG, fraction, fracrepetition, algrepetitions)




# n = 10000
# reads = collect(1:n)
# sets = [Set(sample(reads, rand(1:length(reads)), replace=false)) for _ in 1:n/10]

# jaccardindex(s1, s2) = length(s1 ∩ s2) / length(s1 ∪ s2)

# j_m = Matrix{Float64}(undef, length(sets), length(sets))

# for x in eachindex(sets)
#     s1 = sets[x]
#     for y in eachindex(sets)
#         s2 = sets[y]
#         j_m[x, y] = jaccardindex(s1, s2)
#     end
# end



