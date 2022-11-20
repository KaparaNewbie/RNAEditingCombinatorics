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






"""
    distances(M, AA_groups)

Create a symmaetrical distances matrix `Δ` which measures the distance between any two 
`rowᵢ, rowⱼ ∈ M`.  
The distance between `rowᵢ` to `rowⱼ` is determined by the number of corresponding columns 
in which the two contain no possible amino acids `AAᵦ` and `AAᵧ`, respectively, 
such that both `AAᵦ` and `AAᵧ` share the same classification in `AA_groups`.
"""
function distances(M::Matrix{Set{AminoAcid}}, AA_groups::Dict{AminoAcid,String})
    nrows = size(M, 1)
    Δmax = size(M, 2)  # maximal possible difference between any two proteins
    uint = smallestuint(Δmax)  # type of each cell in the final distances matrix `Δ`
    AAsets = ThreadsX.Set([x for row ∈ eachrow(M) for x ∈ row])
    # sets with no two possible AAs that share a similar classification
    distinctAAsets = ThreadsX.Set(
        [
        (x, y)
        for x ∈ AAsets for y ∈ AAsets
        if !anysimilarity(x, y, AA_groups)
    ]
    )
    # for each `rowᵢ`, calc `δᵢ`which is a vector of distances to each `rowⱼ` where `i <= j`
    Δs = tcollect(i_js_distances(M, distinctAAsets, uint, i) for i ∈ 1:nrows)
    # merge the triangular array `Δs` into a symmaetrical distances matrix `Δ`
    Δ = SymmetricPacked(reduce(vcat, Δs))
    return Δ
end




"""
    distances(M, substitutionmatrix, minsimilarityscore, similarityvalidator)

Create a symmaetrical distances matrix `Δ` which measures the distance between any two 
`rowᵢ, rowⱼ ∈ M`.  
The distance between `rowᵢ` to `rowⱼ` is determined by the number of corresponding columns 
in which the two share no possible amino acid.
"""
function distances(
    M::Matrix{Set{AminoAcid}},
    substitutionmatrix::SubstitutionMatrix{AminoAcid,Int64}, minsimilarityscore::Int64, similarityvalidator::Function
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
        if !anysimilarity(x, y, substitutionmatrix, minsimilarityscore, similarityvalidator)
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

Create a symmaetrical distances matrix `Δ` which measures the distance between any two 
`rowᵢ, rowⱼ ∈ M`.  
The distance between `rowᵢ` to `rowⱼ` is determined by the number of corresponding columns 
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
    M::Matrix{Set{AminoAcid}}, distinctAAsets, uint, i
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
    maxmainthreads, outdir,
    substitutionmatrix::Union{SubstitutionMatrix,Nothing},
    minsimilarityscore::Int64,
    similarityvalidator::Function,
    useAAgroups::Bool
)
    distinctdf = prepare_distinctdf(distinctfile, delim, innerdelim, truestrings, falsestrings)

    allprotsdf, firstcolpos = prepare_allprotsdf!(
        allprotsfile, delim, innerdelim, truestrings, falsestrings, firstcolpos
    )

    # the possible amino acids each protein has in each position
    M = Matrix(allprotsdf[:, firstcolpos:end])
    # the distances between any two proteins according to `M`
    # Δ = distances(M) 
    Δ = begin
        if substitutionmatrix !== nothing
            distances(M, substitutionmatrix, minsimilarityscore, similarityvalidator)
        elseif useAAgroups
            distances(M, AA_groups)
        else
            distances(M)
        end
    end

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
    outfile = joinpath(abspath(outdir), "$samplename.DistinctUniqueProteins.ExpressionLevels$postfix_to_add.csv")
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
        "--distinctfiles"
        help = "One or more csv files representing distinct unique proteins."
        nargs = '+'
        action = :store_arg
        required = true
        "--allprotsfiles"
        help = "corresponding csv files representing unique proteins (w.r.t. `distinctfiles`)."
        nargs = '+'
        action = :store_arg
        required = true
        "--samplenames"
        nargs = '+'
        action = :store_arg
        required = true
        "--postfix_to_add"
        help = "Add `postfix` to output files' names, e.g., `\$sample.DistinctUnique{Reads,Proteins}\$postfix.\$time.csv`."
        default = ""
        "--firstcolpos"
        help = "Int location of the first editing position column of each file in `infiles`. As of now, should be 9 for `Reads` and 15 for `Proteins`."
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

    end
    return parse_args(s)
end



function main(
    distinctfiles, allprotsfiles, samplenames, postfix_to_add,
    firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
    maxmainthreads, outdir,
    gcp, shutdowngcp,
    substitutionmatrix::Union{SubstitutionMatrix,Nothing},
    minsimilarityscore::Int64,
    similarityvalidator::Function,
    useAAgroups::Bool
)
    # run each sample using the `run_sample` function
    for (distinctfile, allprotsfile, samplename) ∈ zip(distinctfiles, allprotsfiles, samplenames)
        run_sample(
            distinctfile, allprotsfile, samplename, postfix_to_add,
            firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
            maxmainthreads, outdir,
            substitutionmatrix, minsimilarityscore, similarityvalidator, useAAgroups
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
    postfix_to_add = parsedargs["postfix_to_add"]

    firstcolpos = parsedargs["firstcolpos"]
    delim = parsedargs["delim"]
    innerdelim = parsedargs["innerdelim"]
    truestrings = parsedargs["truestrings"]
    falsestrings = parsedargs["falsestrings"]

    fractions = parsedargs["fractions"]
    fractions = fractions isa Vector{Float64} ? fractions : parse.(Float64, fractions)

    maxmainthreads = parsedargs["maxmainthreads"]
    outdir = parsedargs["outdir"]

    gcp = parsedargs["gcp"]
    shutdowngcp = parsedargs["shutdowngcp"]

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

    # run
    main(
        distinctfiles, allprotsfiles, samplenames, postfix_to_add,
        firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
        maxmainthreads, outdir,
        gcp, shutdowngcp,
        substitutionmatrix, minsimilarityscore, similarityvalidator, useAAgroups
    )
end


if abspath(PROGRAM_FILE) == @__FILE__
    CLI_main()
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

# samplenames = [
#     "RUSC2_MOUSE",
#     "TRIM2_BOVIN",
#     "CA2D3_MOUSE",
#     "ABL_DROME",
#     "DGLA_HUMAN",
#     "K0513_MOUSE",
#     "KCNAS_DROME",
#     "ACHA4_MOUSE",
#     "ANR17_HUMAN",
#     "TWK7_CAEEL",
#     "SCN1_HETBL",
#     "CACB2_RABIT"
# ]
# chroms = [
#     "comp141881_c0_seq3",
#     "comp141044_c0_seq2",
#     "comp140439_c0_seq1",
#     "comp126362_c0_seq1",
#     "comp141517_c0_seq1",
#     "comp141840_c0_seq2",
#     "comp141640_c0_seq1",
#     "comp140987_c3_seq1",
#     "comp140910_c2_seq1",
#     "comp136058_c0_seq1",
#     "comp141378_c0_seq7",
#     "comp141158_c1_seq2"
# ]
# allprotsfiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.$chrom.unique_proteins.csv"
#     for chrom in chroms
# ]
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina"

# firstcolpos = 15

# delim = "\t"
# innerdelim = ","

# truestrings = ["TRUE", "True", "true"]
# falsestrings = ["FALSE", "False", "false"]

# fractions = [1.0]

# maxmainthreads = 30



# main(
#     distinctfiles, allprotsfiles, samplenames,
#     firstcolpos, delim, innerdelim, truestrings, falsestrings, fractions,
#     maxmainthreads, outdir
# )


