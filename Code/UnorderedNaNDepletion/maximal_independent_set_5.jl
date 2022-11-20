# https://discourse.julialang.org/t/julia-python-equivalent-of-main/35433
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    using Distributed
else
    using BenchmarkTools
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
using BioAlignments
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
    const to = TimerOutput() # Create a TimerOutput, this is the main type that keeps track of everything.
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
        help = "Use this substitution matrix as a stricter criteria for determination of distinct AAs. Use in conjuction with `datatype == Proteins`. Not compatible with `AA_groups`."
        "--minsimilarityscore"
        help = "See `similarityvalidator` below."
        arg_type = Int
        default = 0
        "--similarityvalidator"
        help = "Use this opeartor to determine similarty of AA change, e.g., whether `5 >= minsimilarityscore`."
        arg_type = Symbol
        default = :(>=)
        "--useAAgroups"
        help = "Classify AAs by groups (polar/non-polar/positive/negative) as a stricter criteria for determination of distinct AAs. Use in conjuction with `datatype == Proteins`. Not compatible with `substitutionmatrix`."
        action = :store_true
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

    _substitutionmatrix = parsedargs["substitutionmatrix"]
    if _substitutionmatrix !== nothing # it's a string with a matrix name
        # substitutionmatrix = BioAlignments.load_submat(BioSymbols.AminoAcid, uppercase(_substitutionmatrix))
        substitutionmatrix = BioAlignments.load_submat(BioSequences.AminoAcid, uppercase(_substitutionmatrix))
    else
        substitutionmatrix = nothing
    end

    minsimilarityscore = parsedargs["minsimilarityscore"]
    similarityvalidator = eval(parsedargs["similarityvalidator"])
    useAAgroups = parsedargs["useAAgroups"]





    # infiles = ["D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"]
    # delim = "\t"
    # postfix_to_remove = ".aligned.sorted.MinRQ998.unique_proteins.csv"
    # prefix_to_remove = ""
    # postfix_to_add = ".SmallCloudTest"
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
    samplesnames = map(x -> replace(splitdir(x)[2], prefix_to_remove => "", postfix_to_remove => ""), infiles)

    # logfile = joinpath(outdir, "log.$(writingtime()).txt")
    # run(`touch $logfile`)
    # @everywhere setlogger($logfile) # https://discourse.julialang.org/t/copying-variable-to-all-remote-processes/26518/4?u=kaparanewbie

    for x ∈ eachindex(infiles)  # could be eachindex(samplesnames) as well
        infile = infiles[x]
        samplename = samplesnames[x]
        @timeit to "run_sample" run_sample(
            infile, delim, samplename, idcol, firstcolpos, datatype, outdir, postfix_to_add,
            fracstep, maxfrac, fracrepetitions, algrepetitions, testfraction, randseed,
            run_solve_threaded, sortresults, algs,
            substitutionmatrix, minsimilarityscore, similarityvalidator,
            useAAgroups
        )
    end

    # Print the timings in the default way
    show(to)

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
    minsimilarityscore::Int64,
    similarityvalidator::Function,
    useAAgroups::Bool
)
    @info "$(loggingtime())\trun_sample" infile samplename

    df, firstcolpos = try
        @timeit to "preparedf!" preparedf!(
            infile, delim, datatype, idcol, firstcolpos,
            testfraction, randseed
        )
    catch e
        @warn "$(loggingtime())\tpreparedf! failed for $infile" e
        return
    end


    # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
    G = try
        if substitutionmatrix !== nothing
            @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol, substitutionmatrix, minsimilarityscore, similarityvalidator; firstcolpos)
        elseif useAAgroups
            @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol, AA_groups; firstcolpos)
        else
            @timeit to "indistinguishable_rows" indistinguishable_rows(df, idcol; firstcolpos)
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
    @timeit to "fracrepetitions_inputs" fracrepetitions_inputs = prep_pmap_input(fracstep, maxfrac, df, fracrepetitions)

    @timeit to "run_fracrepetition" fracrepetitions_results = pmap(
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
            @warn "$(loggingtime())\trun_fracrepetition failed" fraction fracrepetition e
            missing
        end
    end

    @timeit to "`fracrepetitions_results` -> sorted `results`" begin
        results::DataFrame = vcat(skipmissing(fracrepetitions_results)...)
        sort!(results, ["Fraction", "FractionRepetition", "Algorithm", "AlgorithmRepetition"])
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


# infile = "D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
# delim = "\t"
# postfix_to_remove = ".aligned.sorted.MinRQ998.unique_proteins.csv"
# prefix_to_remove = ""
# postfix_to_add = ".SmallCloudTest"
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



# df, firstcolpos = preparedf!(
#     infile, delim, datatype, idcol, firstcolpos,
#     testfraction, randseed
# )


# substitutionmatrix = BLOSUM62
# minsimilarityscore = 0
# similarityvalidator = ≥






# # G is the main neighborhood matrix and created only once; samples can create subgraphs induced by it
# G1 = indistinguishable_rows(df, idcol; firstcolpos)
# G2 = indistinguishable_rows(df, idcol, substitutionmatrix, minsimilarityscore, similarityvalidator; firstcolpos)
# G3 = indistinguishable_rows(df, idcol, AA_groups; firstcolpos)

# numedges1 = sum(length.(values(G1)))
# numedges2 = sum(length.(values(G2)))
# numedges3 = sum(length.(values(G3)))

# e = :("BLOSUM45")
# eval(e)
# submat = eval("BLOSUM45")



# BioAlignments.load_submat(BioSymbols.AminoAcid, "BLOSUM45")

# e1 = :(">=")
# e2 = :(>=)

# typeof(e1)
# typeof(e2)

# Symbol(e1)
# stringop = ">="
# e3 = :($stringop)

# eval(e1)
# eval(e2)


# eval(Symbol(e1))(5, 3)

# eval(:($">="))

# eval(:(>=))



# op = 