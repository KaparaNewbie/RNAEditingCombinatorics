# TODO complete wrapping this script to be run from the CLI


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



"""
    distances(M)

Create a symmaetrical distances matrix `Δ` which measures the distance between any two 
`rowᵢ, rowⱼ ∈ M`.  
The distance between any `rowᵢ` to `rowⱼ` is determined by to the number of corresponding columns 
in which the two share no possible amino acid.
"""
function distances(M::Matrix{Set{AminoAcid}})
    nrows = size(M, 1)
    Δmax = size(M, 2)  # maximal possible difference between any two proteins
    uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
    emptyAAintersections = ThreadsX.Set(
        [(x, y)
         for x ∈ AAsets
         for y ∈ AAsets
         if x ∩ y == Set()]
    )
    # for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
    Δs = tcollect(i_js_distances(M, emptyAAintersections, uint, i) for i ∈ 1:nrows)
    # merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
    Δ = SymmetricPacked(reduce(vcat, Δs))
    return Δ
end


"""
    i_js_distances(M, emptyAAintersections, uint, i)

Return `δᵢ` which is a vector of distances between `rowᵢ` to every `rowⱼ` in `M`, where `i <= j <= size(M, 1)`.  
A distance between `rowᵢ` to a `rowⱼ` is determined by to the number of corresponding positions in which the 
two share no possible amino acid (naturally, the distance between `rowᵢ` to itself is 0). 
All empty intersections between any two cells are given at `emptyAAintersections`.
`uint` is the type of unsigned int of the distances.
"""
function i_js_distances(
    M::Matrix{Set{AminoAcid}}, emptyAAintersections, uint, i
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
            # the two rows of proteins have no amino acid they possibly share
            if (x, y) ∈ emptyAAintersections
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
    distinctfile, allprotsfile, samplename,
    firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
    maxmainthreads, outdir
)
    distinctdf = prepare_distinctdf(distinctfile, delim, innerdelim, truestrings, falsestrings)

    allprotsdf, firstcolpos = prepare_allprotsdf!(
        allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos
    )

    M = Matrix(allprotsdf[:, firstcolpos:end]) # the possible amino acids each protein has in each position
    Δ = distances(M) # the distances between any two proteins according to `M`

    # considering only solutions of desired fractions
    solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[:, "Index"]

    # basesize = Int(round(Threads.nthreads()))
    # maxmainthreads = Int(round(Threads.nthreads() / 5))
    # basesize = Int(round(Threads.nthreads() / 5))
    allsubsolutions = collect(Iterators.partition(solutions, maxmainthreads))
    results = tcollect(
        additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
        for subsolutions ∈ allsubsolutions
    )
    finalresults = vcat(Iterators.flatten(results)...)

    # save the results
    outfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels.csv")
    CSV.write(outfile, finalresults; delim)
end


function additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, solutions)
    results = map(solutions) do solution
        one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, solution)
    end
    return results
    # vcat(results)
end


function prepare_distinctdf(distinctfile, delim, innerdelim, truestrings, falsestrings)
    distinctdf = DataFrame(CSV.File(distinctfile; delim, truestrings, falsestrings))
    transform!(distinctdf, :UniqueSamples => (x -> split.(x, innerdelim)) => :UniqueSamples)
    distinctdf[!, "Index"] = collect(1:size(distinctdf, 1))
    return distinctdf
end


toAAset(x, innerdelim) = Set(map(aa -> convert(AminoAcid, only(aa)), split(x, innerdelim)))


function prepare_allprotsdf!(allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos)
    allprotsdf = DataFrame(CSV.File(allprotsfile; delim, truestrings, falsestrings))
    allprotsdf = hcat(allprotsdf[:, begin:firstcolpos-1], toAAset.(allprotsdf[:, firstcolpos:end], innerdelim))

    insertcols!(allprotsdf, firstcolpos, :Index => 1:size(allprotsdf, 1))
    firstcolpos += 1
    insertcols!(allprotsdf, firstcolpos, :AdditionalEqualSupportingReads => 0.0)
    firstcolpos += 1
    insertcols!(allprotsdf, firstcolpos, :AdditionalWeightedSupportingReads => 0.0)
    firstcolpos += 1

    return allprotsdf, firstcolpos
end


# solution = 901

# unchosendf = one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, solution)

# describe(unchosendf[!, ["NumOfReads", "AdditionalEqualSupportingReads", "AdditionalWeightedSupportingReads", "TotalEqualSupportingRead", "TotalWeightedSupportingRead"]])
# sum.(eachcol(unchosendf[!, ["NumOfReads", "AdditionalEqualSupportingReads", "AdditionalWeightedSupportingReads", "TotalEqualSupportingRead", "TotalWeightedSupportingRead"]]))

function one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, solution)

    solutionrow = distinctdf[solution, :]

    prots_in_solution = solutionrow["UniqueSamples"]

    allprotsdf = allprotsdf[:, begin:firstcolpos-1]

    # sum(allprotsdf[:, "NumOfReads"])

    chosendf = filter("Protein" => x -> x ∈ prots_in_solution, allprotsdf)
    unchosendf = filter("Protein" => x -> x ∉ prots_in_solution, allprotsdf)

    insertcols!(chosendf, 3, "#Solution" => solutionrow["Index"])
    insertcols!(chosendf, 4, "Fraction" => solutionrow["Fraction"])
    insertcols!(chosendf, 5, "FractionRepetition" => solutionrow["FractionRepetition"])
    insertcols!(chosendf, 6, "Algorithm" => solutionrow["Algorithm"])
    insertcols!(chosendf, 7, "AlgorithmRepetition" => solutionrow["AlgorithmRepetition"])

    chosenindices = chosendf[!, "Index"] # indices of chosen proteins in the complete Δ matrix

    for unchosenprot ∈ eachrow(unchosendf)  # a row of unchosen protein relative to the chosen proteins in the solution

        # unchosenprot = unchosendf[1000, :]  # an example row of unchosen protein

        unchosenprot_index = unchosenprot["Index"] # the row number of the unchosen protein in the distances matrix Δ

        unchosenprot_distances = Δ[unchosenprot_index, chosenindices] # the distances between the unchosen protein `unchosenprot` to all chosen proteins

        unchosenprot_distances_argmins = allargmins(unchosenprot_distances) # indices of minimum distances

        minchosendf = chosendf[unchosenprot_distances_argmins, :] # chosen proteins with minimum distance to the unchosen protein `unchosenprot`

        unchosenreads = unchosenprot["NumOfReads"] # the number of reads `unchosenprot` divides between the chosen proteins that are closest to it

        equal_addition = unchosenreads / size(minchosendf, 1)
        weighted_addition = unchosenreads .* minchosendf[!, "NumOfReads"] ./ sum(minchosendf[!, "NumOfReads"])

        # sum(chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"])
        # sum(chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"])

        chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"] .+= equal_addition
        chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"] .+= weighted_addition

        # sum(chosendf[unchosenprot_distances_argmins, "AdditionalEqualSupportingReads"])
        # sum(chosendf[unchosenprot_distances_argmins, "AdditionalWeightedSupportingReads"])

    end

    chosendf[!, "TotalEqualSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalEqualSupportingReads"]
    chosendf[!, "TotalWeightedSupportingReads"] .= chosendf[!, "NumOfReads"] .+ chosendf[!, "AdditionalWeightedSupportingReads"]

    return chosendf

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



function main(
    distinctfiles, allprotsfiles, samplenames,
    firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
    maxmainthreads, outdir
)

    # todo (1) create a CLI using Comonicon.jl https://github.com/comonicon/Comonicon.jl 

    # (2) run each sample using the `run_sample` function
    for (distinctfile, allprotsfile, samplename) ∈ zip(distinctfiles, allprotsfiles, samplenames)
        run_sample(
            distinctfile, allprotsfile, samplename,
            firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
            maxmainthreads, outdir
        )
    end

end



# distinctfiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueProteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueProteins.csv"
# ]

# allprotsfiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"
# ]

# samplenames = ["GRIA", "PCLO"]
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2"

distinctfiles = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141881_c0_seq3.DistinctUniqueProteins.12.07.2022-20:54:38.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141044_c0_seq2.DistinctUniqueProteins.13.07.2022-06:33:23.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140439_c0_seq1.DistinctUniqueProteins.12.07.2022-22:51:22.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp126362_c0_seq1.DistinctUniqueProteins.15.07.2022-06:11:18.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141517_c0_seq1.DistinctUniqueProteins.14.07.2022-07:43:15.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141840_c0_seq2.DistinctUniqueProteins.13.07.2022-20:30:25.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141640_c0_seq1.DistinctUniqueProteins.12.07.2022-19:44:02.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140987_c3_seq1.DistinctUniqueProteins.18.07.2022-07:50:43.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp140910_c2_seq1.DistinctUniqueProteins.13.07.2022-16:15:35.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp136058_c0_seq1.DistinctUniqueProteins.21.07.2022-07:57:53.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141378_c0_seq7.DistinctUniqueProteins.19.07.2022-08:12:24.csv",
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/comp141158_c1_seq2.DistinctUniqueProteins.13.07.2022-01:54:59.csv",
]

samplenames = [
    "RUSC2_MOUSE",
    "TRIM2_BOVIN",
    "CA2D3_MOUSE",
    "ABL_DROME",
    "DGLA_HUMAN",
    "K0513_MOUSE",
    "KCNAS_DROME",
    "ACHA4_MOUSE",
    "ANR17_HUMAN",
    "TWK7_CAEEL",
    "SCN1_HETBL",
    "CACB2_RABIT"
]
chroms = [
    "comp141881_c0_seq3",
    "comp141044_c0_seq2",
    "comp140439_c0_seq1",
    "comp126362_c0_seq1",
    "comp141517_c0_seq1",
    "comp141840_c0_seq2",
    "comp141640_c0_seq1",
    "comp140987_c3_seq1",
    "comp140910_c2_seq1",
    "comp136058_c0_seq1",
    "comp141378_c0_seq7",
    "comp141158_c1_seq2"
]
allprotsfiles = [
    "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.$chrom.unique_proteins.csv"
    for chrom in chroms
]
outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina"

firstcolpos = 15

delim = "\t"
innerdelim = ","

truestrings = ["TRUE", "True", "true"]
falsestrings = ["FALSE", "False", "false"]

fractions = [1.0]

maxmainthreads = 30



main(
    distinctfiles, allprotsfiles, samplenames,
    firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
    maxmainthreads, outdir
)


# (1) this is what sent to `additional_assignments` by `main`

# distinctfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueProteins.csv"
# allprotsfile = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv"

# firstcolpos = 15

# delim = "\t"
# innerdelim = ","

# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]

# fractions = [1.0]

# maxmainthreads = Int(round(Threads.nthreads() / 5))


# # (2) this is run by `additional_assignments`

# distinctdf = prepare_distinctdf(distinctfile, delim, innerdelim, truestrings, falsestrings)

# allprotsdf, firstcolpos = prepare_allprotsdf!(
#     allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos
# )

# M = Matrix(allprotsdf[:, firstcolpos:end]) # the possible amino acids each protein has in each position
# Δ = distances(M) # the distances between any two proteins according to `M`



# # considering only solutions of desired fractions
# # solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[:, "Index"]
# # basesize = Int(round(Threads.nthreads() / 5))



# # ascsolutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[1:10:99, :]
# # descsolutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[6:10:99, :]
# # ascsolutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[1:10:99, "Index"]
# # descsolutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[6:10:99, "Index"]

# # solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[1:4, "Index"]
# # solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[1:10:99, "Index"]
# solutions = subset(distinctdf, "Fraction" => x -> x .∈ fractions)[collect(1:10:99) ∪ collect(6:10:99), "Index"]
# basesize = 1

# allsubsolutions = collect(Iterators.partition(solutions, basesize))

# results = tcollect(
#     additional_assignments(distinctdf, allprotsdf, firstcolpos, Δ, subsolutions)
#     for subsolutions ∈ allsubsolutions
# )

# finalresults = vcat(Iterators.flatten(results)...)
# # finalresults = vcat(results)


# samplename = "GRIA"
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2"
# outfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels.csv")
# CSV.write(outfile, finalresults; delim)






# r1 = one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, 901)
# r2 = one_solution_additional_assignment(distinctdf, allprotsdf, firstcolpos, Δ, 905) 

# r1[!, "#Solution"] .= string.(r1[!, "#Solution"])
# r2[!, "#Solution"] .= string.(r2[!, "#Solution"])

# r = vcat(r1, r2)


# r[!, "#Solution"] .= string.(r[!, "#Solution"])

# using AlgebraOfGraphics
# using CairoMakie



# plt = data(r1) * mapping(:TotalEqualSupportingRead) * histogram(bins=20)
# fg = draw(plt)


# plt4 = data(r1) * mapping(:TotalEqualSupportingRead) * histogram(bins=20)
# fg = draw(plt)


# plt2 = data(r) * mapping(:TotalEqualSupportingRead, col="#Solution") * histogram(bins=20)
# fg2 = draw(plt2)


# plt3 = data(r) * mapping(:TotalEqualSupportingRead, col="#Solution") * histogram(bins=20)
# fg23 = draw(plt23)


# ###

# resolution = (800, 600)
# ff = Figure(; resolution)
# ax1 = Axis(ff[1, 1])
# ax2 = Axis(ff[1, 2])

# plt1 = data(r) * mapping(:TotalEqualSupportingRead, color="#Solution", stack="#Solution") * histogram(bins=20)
# plt2 = data(r) * mapping(:TotalWeightedSupportingRead, color="#Solution", stack="#Solution") * histogram(bins=20)

# grid = draw!(ax1, plt1)
# grid = draw!(ax2, plt2)

# legend!(ff[1, 3], grid)

# ff



# resolution = (800, 600)
# ff = Figure(; resolution)

# ax1 = Axis(ff[1, 1])
# ax2 = Axis(ff[1, 2])
# ax3 = Axis(ff[2, 1])
# ax4 = Axis(ff[2, 2])

# plt1 = data(r1) * mapping(:TotalEqualSupportingRead) * histogram(bins=20)
# plt2 = data(r1) * mapping(:TotalWeightedSupportingRead) * histogram(bins=20)
# plt3 = data(r2) * mapping(:TotalEqualSupportingRead) * histogram(bins=20)
# plt4 = data(r2) * mapping(:TotalWeightedSupportingRead) * histogram(bins=20)

# grid = draw!(ax1, plt1)
# grid = draw!(ax2, plt2)
# grid = draw!(ax3, plt3)
# grid = draw!(ax4, plt4)

# # legend!(ff[1, 3], grid)

# ff




# using AlgebraOfGraphics: density

# resolution = (800, 600)
# ff = Figure(; resolution)
# ax1 = Axis(ff[1, 1])
# ax2 = Axis(ff[1, 2])

# plt1 = data(r) * mapping(:TotalEqualSupportingRead, color="#Solution") * density(bandwidth=0.5)
# plt2 = data(r) * mapping(:TotalWeightedSupportingRead, color="#Solution") * density(bandwidth=0.5)

# grid = draw!(ax1, plt1)
# grid = draw!(ax2, plt2)

# legend!(ff[1, 3], grid)

# ff
###

