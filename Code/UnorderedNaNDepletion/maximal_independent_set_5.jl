# https://discourse.julialang.org/t/julia-python-equivalent-of-main/35433
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    using Distributed
else
    using BenchmarkTools
    using Distributed
end


n_workers = length(workers())
if n_workers > 1
    # remove previous workers
    rmprocs(workers())
    # add new one, but with limited number of threads
    addprocs(
        # n_workers  = 2
        n_workers,
        exeflags=[
            "--threads=$(Int(round(Threads.nthreads() / n_workers)))",
            "--project"
        ]
    )
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
using TimerOutputs
using BioAlignments  # right now it's a local branch from my own repo, 
using BioSymbols


include(joinpath(@__DIR__, "consts.jl")) # for ∅ & AA_groups
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
    using TimerOutputs
    # using BioSymbols
    using BioSequences
    # const to = TimerOutput() # Create a TimerOutput, this is the main type that keeps track of everything.
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
        help = "Add `postfix` to output files' names, e.g., `\$sample.DistinctUnique{Reads,Proteins}\$postfix.\$time.csv`."
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

    substitutionmatrix = eval(parsedargs["substitutionmatrix"])
    similarityscorecutoff = parsedargs["similarityscorecutoff"]
    similarityvalidator = eval(parsedargs["similarityvalidator"])
    aagroups = eval(parsedargs["aagroups"])


    # infiles = ["D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"]
    # delim = "\t"
    # postfix_to_remove = ".aligned.sorted.MinRQ998.unique_proteins.csv"
    # prefix_to_remove = ""
    # postfix_to_add = ".GRANTHAM1974-100"
    # idcol = "Protein"
    # firstcolpos = 15
    # datatype = "Proteins"
    # outdir = "D.pealeii/MpileupAndTranscripts/RQ998.2"
    # fracstep = 0.5
    # maxfrac = 1.0
    # fracrepetitions = 5
    # algrepetitions = 2
    # testfraction = 0.01
    # randseed = 1892
    # run_solve_threaded = false
    # sortresults = false
    # algs = ["Ascending", "Descending"]
    # gcp = false
    # shutdowngcp = false


    @info "$(loggingtime())\tmain" delim prefix_to_remove postfix_to_remove postfix_to_add idcol firstcolpos datatype outdir fracstep maxfrac fracrepetitions algrepetitions testfraction randseed run_solve_threaded sortresults algs gcp shutdowngcp

    # do the thing for each file, in ascending order of file size

    sort!(infiles, by=(x -> filesize(x)))
    if prefix_to_remove != "" && postfix_to_remove != ""
        samplesnames = map(x -> replace(splitdir(x)[2], prefix_to_remove => "", postfix_to_remove => ""), infiles)
    elseif prefix_to_remove != ""
        samplesnames = map(x -> replace(splitdir(x)[2], prefix_to_remove => ""), infiles)
    elseif postfix_to_remove != ""
        samplesnames = map(x -> replace(splitdir(x)[2], postfix_to_remove => ""), infiles)
    else
        samplesnames = map(x -> splitdir(x)[2], infiles)
    end

    # logfile = joinpath(outdir, "log.$(writingtime()).txt")
    # run(`touch $logfile`)
    # @everywhere setlogger($logfile) # https://discourse.julialang.org/t/copying-variable-to-all-remote-processes/26518/4?u=kaparanewbie

    for x ∈ eachindex(infiles)  # could be eachindex(samplesnames) as well
        infile = infiles[x]
        samplename = samplesnames[x]
        # @timeit to "run_sample" run_sample(
        #     infile, delim, samplename, idcol, firstcolpos, datatype, outdir, postfix_to_add,
        #     fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed,
        #     run_solve_threaded, sortresults, algs,
        #     substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups
        # )
        run_sample(
            infile, delim, samplename, idcol, firstcolpos, datatype, outdir, postfix_to_add,
            fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed,
            run_solve_threaded, sortresults, algs,
            substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups
        )
    end

    # # Print the timings in the default way
    # show(to)

    # shutdown gcp vm
    gcp && shutdowngcp && run(`sudo shutdown`) # https://cloud.google.com/compute/docs/shutdownscript
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
    algs::Vector{String},
    substitutionmatrix::Union{SubstitutionMatrix,Nothing},
    similarityscorecutoff::Int64,
    similarityvalidator::Function,
    aagroups::Union{Dict{AminoAcid,String},Nothing}
)
    @info "$(loggingtime())\trun_sample" infile samplename

    df, firstcolpos = try
        # @timeit to "preparedf!" preparedf!(
        #     infile, delim, datatype, idcol, firstcolpos,
        #     testfraction, randseed
        # )
        preparedf!(
            infile, delim, datatype, idcol, firstcolpos,
            testfraction, randseed
        )
    catch e
        @warn "$(loggingtime())\tpreparedf! failed for $infile" e
        return
    end

    if length(size(df)) > 1  # df has 2 dimensions, i.e., not a row of a single protein

        # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
        G = try
            if substitutionmatrix !== nothing
                # @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol, substitutionmatrix, similarityscorecutoff, similarityvalidator; firstcolpos)
                indistinguishable_rows(df, idcol, substitutionmatrix, similarityscorecutoff, similarityvalidator; firstcolpos)
            elseif aagroups !== nothing
                # @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol, aagroups; firstcolpos)
                indistinguishable_rows(df, idcol, aagroups; firstcolpos)
            else
                # @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol; firstcolpos)
                indistinguishable_rows(df, idcol; firstcolpos)
            end
        catch e
            @warn "$(loggingtime())\tindistinguishable_rows failed for $infile" e
            return
        end
        ArrG = @DArray [G for _ ∈ 1:1]  # move G across processes on a distributed array in order to save memory

        # having built G, we only need to keep the reads and unique reads/proteins they support
        select!(df, idcol)
        GC.gc() # free up memory, just in case


        # define input parameters for pmap
        # @timeit to "fracrepetitions_inputs" fracrepetitions_inputs = prep_pmap_input(fracstep, maxfrac, df, fracrepetitions)
        fracrepetitions_inputs = prep_pmap_input(fracstep, maxfrac, df, fracrepetitions)

        # @timeit to "run_fracrepetition" fracrepetitions_results = pmap(
        fracrepetitions_results = pmap(
            fracrepetitions_inputs;
            retry_delays=ExponentialBackOff(n=3, first_delay=5, max_delay=1000)
        ) do (fraction, nsamplerows, fracrepetition)
            try
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
            catch e
                @warn "$(loggingtime())\trun_fracrepetition failed" infile fraction fracrepetition e
                missing
            end
        end

        # @timeit to "`fracrepetitions_results` -> sorted `results`" begin
        #     results::DataFrame = vcat(skipmissing(fracrepetitions_results)...)
        #     sort!(results, ["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
        # end
        results::DataFrame = vcat(skipmissing(fracrepetitions_results)...)
        sort!(results, ["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])


    else

        try
            results = emptyresults()
            for alg in algs
                push!(results, [1.0, 1, alg, 1, 1, df[idcol]]; promote=true)
            end
        catch e
            @warn samplename e
            missing
        end
    end

    # write the results to a csv file
    # outfile = abspath(outdir) * "/" * samplename * ".DistinctUnique$datatype.$(writingtime()).csv"
    outfile = joinpath(abspath(outdir), "$samplename.DistinctUnique$datatype$postfix_to_add.$(writingtime()).csv")
    CSV.write(outfile, results; delim)
end


# reset_timer!(to);


main();

# show(to)
# TimerOutputs.complement!(to);
# show(to)





# @benchmark main()


# infile = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq/ProteinsFiles/comp144504_c0_seq1.unique_proteins.csv"
# firstcolpos = 16
# delim = "\t"
# idcol = "Protein"
# datatype = "Proteins"
# testfraction = 1.0
# randseed = 1892

# samplename = "comp179788_c0_seq1"
# outdir = "O.vulgaris/MpileupAndTranscripts/PRJNA791920/IsoSeq"
# postfix_to_add = ""
# fracstep = 0.2
# maxfrac = 1.0
# fracrepetitions = 5
# algrepetitions = 2
# run_solve_threaded = true
# sortresults = false
# algs = ["Ascending", "Descending"]
# substitutionmatrix = nothing
# aagroups = nothing
# similarityscorecutoff = 0
# similarityvalidator = >=


# df, firstcolpos = preparedf!(
#     infile, delim, datatype, idcol, firstcolpos,
#     testfraction, randseed
# )

# G = indistinguishable_rows(df, idcol; firstcolpos)


# ArrG = @DArray [G for _ ∈ 1:1]  # move G across processes on a distributed array in order to save memory

# # having built G, we only need to keep the reads and unique reads/proteins they support
# select!(df, idcol)


# run_sample(
#     infile, delim, samplename, idcol, firstcolpos, datatype, outdir, postfix_to_add,
#     fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed,
#     run_solve_threaded, sortresults, algs,
#     substitutionmatrix, similarityscorecutoff, similarityvalidator, aagroups
# )






# infile = "D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.gz"
# firstcolpos = 15

# delim = "\t"
# # postfix_to_remove = ".aligned.sorted.MinRQ998.unique_proteins.csv"
# # prefix_to_remove = ""
# # postfix_to_add = ".GRANTHAM1974-100"
# idcol = "Protein"

# datatype = "Proteins"
# # outdir = "D.pealeii/MpileupAndTranscripts/RQ998.2"
# # fracstep = 0.5
# # maxfrac = 1.0
# # fracrepetitions = 5
# # algrepetitions = 2
# testfraction = 1.0
# randseed = 1892
# # run_solve_threaded = false
# # sortresults = false
# # algs = ["Ascending", "Descending"]
# # gcp = false
# # shutdowngcp = false



# df, firstcolpos = preparedf!(
#     infile, delim, datatype, idcol, firstcolpos,
#     testfraction, randseed
# )

# G = indistinguishable_rows(df, idcol; firstcolpos)

# describe(df)

# substitutionmatrix = GRANTHAM1974
# similarityscorecutoff = 100
# similarityvalidator = <

# aagroups = AA_groups_Miyata1979


# # # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
# G1 = indistinguishable_rows(df, idcol; firstcolpos)
# G2 = indistinguishable_rows(df, idcol, substitutionmatrix, similarityscorecutoff, similarityvalidator; firstcolpos)
# G3 = indistinguishable_rows(df, idcol, aagroups; firstcolpos)

# S1 = sum(length.(values(G1)))
# S2 = sum(length.(values(G2)))
# S3 = sum(length.(values(G3)))



# results2 = solve(
#     G2,
#     1.0,
#     1,
#     2,
#     true,
#     false,
#     ["Ascending", "Descending"]
# )


# distinctfile1 = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueProteins.csv"
# distinctfile2 = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.GRANTHAM1974-100.05.12.2022-13:58:52.csv"
# distinctfile3 = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.AAgroupsMiyata1979.04.12.2022-18:09:53.csv"

# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]

# innerdelim = ","

# function prepare_distinctdf(distinctfile, delim, innerdelim, truestrings, falsestrings)
#     distinctdf = DataFrame(CSV.File(distinctfile; delim, truestrings, falsestrings))
#     transform!(distinctdf, :UniqueSamples => (x -> split.(x, innerdelim)) => :UniqueSamples)
#     distinctdf[!, "Index"] = collect(1:size(distinctdf, 1))
#     return distinctdf
# end

# distinctdf1 = prepare_distinctdf(distinctfile1, delim, innerdelim, truestrings, falsestrings)
# distinctdf2 = prepare_distinctdf(distinctfile2, delim, innerdelim, truestrings, falsestrings)
# distinctdf3 = prepare_distinctdf(distinctfile3, delim, innerdelim, truestrings, falsestrings)



# distinct1 = distinctdf1[argmax(distinctdf1[:, "NumUniqueSamples"]), "UniqueSamples"]
# distinct2 = distinctdf2[argmax(distinctdf2[:, "NumUniqueSamples"]), "UniqueSamples"]
# distinct3 = distinctdf3[argmax(distinctdf3[:, "NumUniqueSamples"]), "UniqueSamples"]


# function verify(G, distinctsamples)
#     badcouples = []
#     for (i, x) in enumerate(distinctsamples)
#         for y in distinctsamples[i+1:end]
#             if x in G[y] || y in G[x]
#                 append!(badcouples, (x, y))
#             end
#         end
#     end
#     badcouples
# end


# badcouples1 = verify(G1, distinct1)
# badcouples2 = verify(G2, distinct2)
# badcouples3 = verify(G3, distinct3)






# using BioAlignments
# using BioSymbols
# using UnicodePlots

# sm = GRANTHAM1974
# AAs = [
#     AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, 
#     AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V
# ]

# scores = Int64[]
# for (i, x) in enumerate(AAs)
#     for y in AAs[i:end]
#         append!(scores, sm[x, y])
#     end
# end