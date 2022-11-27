# TODO incorperate this calculation of `jaccardmatrix` at the end of `maximal_independent_set_5.jl`
# (or convert to a CLI app, altough it's less preferable IMO)


using CSV
using DataFrames
using StatsBase  # for StatsBase.sample
import Base.Threads.@spawn
using BenchmarkTools



# n = 1000
# reads = collect(1:n)
# sets = [Set(sample(reads, rand(1:length(reads)), replace=false)) for _ in 1:n/10]


# n = 5
# reads = collect(1:n)
# sets = [Set(sample(reads, rand(1:length(reads)), replace=false)) for _ in 1:n]


jaccardindex(s1, s2) = length(s1 ∩ s2) / length(s1 ∪ s2)


function jaccardmatrix(sets)
    # intiate a matrix with undefined values
    M = Matrix{Float64}(undef, length(sets), length(sets))
    # calculate jaccard index of each two sets and save the result 
    # in the appropriate cell
    Threads.@threads for x ∈ eachindex(sets)
        s1 = sets[x]
        for y ∈ eachindex(sets)
            s2 = sets[y]
            M[x, y] = jaccardindex(s1, s2)
        end
    end
    # return finished matrix
    return M
end



function main()
    
end



# infiles = [
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA.AllRows.DistinctUniqueProteins.csv",
#     "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/PCLO.AllRows.DistinctUniqueProteins.csv"
# ]

datadir = "/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina"

infiles = [f for f in readdir(datadir, join=true) if occursin("DistinctUnique", f)]

dfs = DataFrame.(CSV.File.(infiles, delim="\t"))

allsets = map(df -> Set.(split.(df[:, "UniqueSamples"], ",")), dfs)

jaccardmatrices = jaccardmatrix.(allsets)

outfiles = map(infiles) do infile
    replace(infile, "DistinctUnique" => "JaccardMatrix")
end


# df = DataFrame(jaccardmatrices[1], :auto)


outdelim = "\t"

for (jaccardmatrix, outfile) ∈ zip(jaccardmatrices, outfiles)
    outdf = DataFrame(jaccardmatrix, :auto)
    rename!(outdf, string.(1:size(outdf, 2)))
    CSV.write(outfile, outdf, delim=outdelim)
end




f(a, b) = a + b
f.([1, 2], [3, 4])
f.([1, 2], 5)