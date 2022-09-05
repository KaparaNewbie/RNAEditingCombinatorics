if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
end

using ArgParse
using CSV
using DelimitedFiles
using DataFrames
import Base.Threads.@spawn
using FASTX


function fastatodict(genome)
    reader = open(FASTA.Reader, genome)
    genebychrom = Dict()
    genepatterns = [r"GN=[^ ]+", r"sp[|][^|]+[|][^|]+"]
    splitsby = ["=", "|"]
    splitsindices = [2, 3]
    for record ∈ reader
        chrom = identifier(record)
        desc = description(record)
        for (pattern, splitby, splitindex) ∈ zip(genepatterns, splitsby, splitsindices)
            gene = match(pattern, desc)
            if gene !== nothing
                gene = split(gene.match, splitby)[splitindex]
                break
            end
        end
        if gene === nothing
            gene = missing
        end
        push!(genebychrom, chrom => gene)
    end
    close(reader)
    return genebychrom
end

function parsecmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--indir"
        help = "A folder withg one or more csv files of short-read matrices created by Ruti."
        required = true
        "--outdir"
        help = "A folder to write the reads and unique reads tables to."
        required = true
        "--genome"
        help = "A fasta file with the genome (actually transcriptome) create by trinity."
        required = true
    end
    return parse_args(s)
end

function main()
    parsedargs = parsecmd()
    indir = abspath(parsedargs["indir"])
    outdir = abspath(parsedargs["outdir"])
    genome = abspath(parsedargs["genome"])
    chrompattern = r"comp[0-9]+_c[0-9]+_seq[0-9]+"
    genebychrom = fastatodict(genome)
    infiles = readdir(indir, join=true)
    sort!(infiles, by=filesize)
    for infile ∈ infiles
        # define outfiles' paths and, if they already exist, skip to the next infile
        samplename = split(basename(infile), ".")[1]
        outfiles = [outdir * "/" * samplename * dftype * ".csv" for dftype ∈ [".reads", ".transcripts"]]
        all(isfile.(outfiles)) && continue
        # make the reads and transcripts (unique reads) data frames
        m = CSV.File(infile; delim="\t") |> Tables.matrix
        # change ruti's notation to mine:
        # missing: 0 => -1
        # no editing: A => 0
        # editing: G => 1
        # sequencing error? C, T => -1 (missing as well)  # todo actually, I didn't consider these cases in my script
        @views replace!(m[:, 2:end], 0 => -1, "0" => -1)
        @views replace!(m[:, 2:end], "G" => 1, "A" => 0, "C" => -1, "T" => -1)
        readsnames = DataFrame(Read=String.(m[:, 1])) # df with reads' names
        # todo check that the serial positions here really match all the position in the genome
        # readsdetails = DataFrame(Int8.(m[:, 2:end]), ["Pos$i" for i ∈ 1:size(m, 2)-1]) # df with reads' details
        # readsdetails = Int8.(m[:, 2:end])
        # readsdetails = DataFrame(readsdetails, ["Pos$i" for i ∈ 1:size(m, 2)-1]) # df with reads' details
        readsdetails = DataFrame(m[:, 2:end], ["Pos$i" for i ∈ 1:size(m, 2)-1]) # df with reads' details
        readsdf = hcat(readsnames, readsdetails)
        # filter uninformative cols
        informativecols = map(col -> 1 ∈ col, eachcol(readsdf))
        informativecols[1] = true # the first column is always informative as it stores the read name
        readsdf = readsdf[:, informativecols]
        # filter uninformative rows
        filter!(row -> 1 ∈ row, readsdf)
        # group by reads data (-1, 0, 1 in each position)
        gdf = DataFrames.groupby(readsdf, Not(:Read))
        transcriptsdf = combine(gdf) do sdf
            DataFrame(NumOfReads=size(sdf, 1))
        end
        groupedreads = [join(sdf.Read, ",") for sdf ∈ gdf]
        insertcols!(transcriptsdf, 1, :Reads => groupedreads)
        select!(transcriptsdf, :Reads, :NumOfReads, Not([:Reads, :NumOfReads]))
        insertcols!(transcriptsdf, 1, :Transcript => 1:size(transcriptsdf, 1))
        # extract gene name by chromosome (the chromosome is present in the name of the file)
        chrom = match(chrompattern, infile).match
        gene = genebychrom[chrom]
        insertcols!(transcriptsdf, 1, :Gene => gene)
        # write dataframes to output
        dfs = [readsdf, transcriptsdf]
        # samplename = split(basename(file), ".")[1]
        # outfiles = [outdir * "/" * samplename * dftype * ".csv" for dftype ∈ [".reads", ".transcripts"]]
        for (df, outfile) ∈ zip(dfs, outfiles)
            CSV.write(outfile, df)
        end
    end
end

main()


###############################################################################

# readsdf = DataFrame(
#     Read = String["r1", "r2", "r3", "r4", "r5", "r6"],
#     Pos1 = Int8[-1, 0, 0, 1, -1, 0],
#     Pos2 = Int8[-1, 0, 0, 1, -1, 0],
#     Pos3 = Int8[-1, 1, 1, 0, -1, 1],
#     Pos4 = Int8[-1, -1, -1, -1, -1, -1],
# )

# # mapcols!(x -> replace(x, -1 => missing), readsdf)
# # allmisscol(col) = Set(col) == Set([-1])

# informativecols = map(col -> 1 ∈ col, eachcol(readsdf))
# informativecols[1] = true # the first column is always informative as it stores the read name
# readsdf = readsdf[:, informativecols]

# filter!(row -> 1 ∈ row, readsdf)


# gdf = groupby(readsdf, "Pos1") 

# for sdf ∈ gdf
#     println(sdf)
# end

# gdf = DataFrames.groupby(df, Not(:Read))
# cdf = combine(gdf) do sdf
#     DataFrame(NumOfReads = size(sdf, 1))
# end
# reads = [join(sdf.Read, ",") for sdf in gdf]   
# insertcols!(cdf, 1, :Reads => reads)
# select!(cdf, :Reads, :NumOfReads, Not([:Reads, :NumOfReads]))
# insertcols!(cdf, 1, :Transcript => 1:size(cdf, 1))

###############################################################################

# genome = "/private7/projects/Combinatorics/D.pealeii/Annotations/December2017/orfs_squ.fa"
# file = "/private7/projects/Combinatorics/D.pealeii/Data/RutisMatrices/unique_base_quality_filtered_MATRIX_2_comp141044_c0_seq2_only_paired_full_comp_from_A.txt.gz"
# outdir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Ruti"
# chrompattern = r"comp[0-9]+_c[0-9]+_seq[0-9]+"
# genebychrom = fastatodict(genome)
# # make the reads and transcripts (unique reads) data frames
# # m = CSV.File(file; buffer_in_memory=true, delim="\t") |> Tables.matrix
# m = CSV.File(file; delim="\t") |> Tables.matrix

# @views replace!(m[:, 2:end], 0 => -1, "0" => -1)
# @views replace!(m[:, 2:end], "G" => 1, "A" => 0, "C" => -1, "T" => -1)

# readsnames = DataFrame(Read = String.(m[:, 1])) # df with reads' names
# readsdetails = DataFrame(Int8.(m[:, 2:end]), ["Pos$i" for i ∈ 1:size(m, 2)-1]) # df with reads' details
# readsdf = hcat(readsnames, readsdetails)

# # filter uninformative cols
# informativecols = map(col -> 1 ∈ col, eachcol(readsdf))
# informativecols[1] = true # the first column is always informative as it stores the read name
# readsdf = readsdf[:, informativecols]

# # filter uninformative rows
# filter!(row -> 1 ∈ row, readsdf)


# gdf = groupby(readsdf, Not(:Read))

# transcriptsdf = combine(gdf) do sdf
#     DataFrame(NumOfReads=size(sdf, 1))
# end

# groupedreads = [join(sdf.Read, ",") for sdf ∈ gdf]
# insertcols!(transcriptsdf, 1, :Reads => groupedreads)
# select!(transcriptsdf, :Reads, :NumOfReads, Not([:Reads, :NumOfReads]))
# insertcols!(transcriptsdf, 1, :Transcript => 1:size(transcriptsdf, 1))
# # extract gene name by chromosome (the chromosome is present in the name of the file)
# chrom = match(chrompattern, file).match
# gene = genebychrom[chrom]
# insertcols!(transcriptsdf, 1, :Gene => gene)
# # write dataframes to output
# dfs = [readsdf, transcriptsdf]
# samplename = split(basename(file), ".")[1]
# outfiles = [outdir * "/" * samplename * dftype * ".csv" for dftype ∈ [".reads", ".transcripts"]]
# for (df, outfile) ∈ zip(dfs, outfiles)
#     CSV.write(outfile, df)
# end

# m2 = CSV.File(f2; buffer_in_memory=true, delim="\t") |> Tables.matrix

# ncols = size(m2, 2)
# dfnames = ["Read"]
# for i in 1:ncols-1
#     push!(dfnames, "Pos$i")
# end

# df2 = DataFrame(m2, dfnames)

# gdf2 = DataFrames.groupby(df2, Not(:Read))

# cdf2 = combine(gdf2) do sdf
#     DataFrame(NumOfReads = size(sdf, 1))
# end

# reads2 = [join(sdf.Read, ",") for sdf in gdf2]   
# insertcols!(cdf2, 1, :Reads => reads2)

# select!(cdf2, :Reads, :NumOfReads, Not([:Reads, :NumOfReads]))

# insertcols!(cdf2, 1, :Transcript => 1:size(cdf2, 1))



# using FASTX

# chrompattern = r"comp[0-9]+_c[0-9]+_seq[0-9]+"
# chrom = match(chrompattern, f2).match

# genome = "/private7/projects/Combinatorics/D.pealeii/Annotations/December2017/orfs_squ.fa"
# reader = open(FASTA.Reader, genome)
# records = [record for record ∈ reader]

# close(reader)

# rec = records[1]

# genepattern = r"GN=[^ ]+"
# gene = match(genepattern, description(rec)).match

# genebychrom = Dict(identifier(rec) => match(genepattern, description(rec)).match for rec in records)

# genebychrom["comp24855_c0_seq1"]


# gene = "blah blah"
# insertcols!(cdf2, 1, :Gene => gene)


# indir = "/private7/projects/Combinatorics/D.pealeii/Data/RutisMatrices"
# infiles = readdir(indir, join=true)
# sort!(infiles, by=filesize)
# for f ∈ infiles
#     println("$(basename(f)), $(Base.format_bytes(filesize(f)))")
# end

