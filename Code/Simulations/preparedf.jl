# using BioSequences # for BioSequences.toAAset
# using DataFrames
# using CSV
# using Random # for Rand.seed
# using StatsBase  # for StatsBase.sample





toAAset(x) = Set(map(aa -> convert(AminoAcid, only(aa)), split(x, ",")))


"""
    preparedf!(infile, delim, removeuneditedrows, allrowsedited, editingcol, datatype, idcol, firstcolpos, testfraction::Union{Float64,Nothing}=nothing, randseed=1892)

Parse `infile` into a DataFrame to be later used for building a graph of (distinct) unique rows (reads / proteins).

# Arguments
- `infile::String`: a csv input file.
- `delim::String`: the delimiter of `infile`.
- `datatype::String`: the data type, either `Reads` or `Proteins`. 
- `idcol::String`: the column name of the column containing the unique row identifier.
- `firstcolpos::Int`: the position of the first column in `infile` containing data about editing in nucleotides / AAs.
- `testfraction::Float64=1.0`: for testing purposes, use only `testfraction < 1.0` of sampled rows. Default `testfraction == 1.0` means that all rows are used.
- `randseed=1892`: the seed for the random number generator in case `testfraction < 1.0`.
"""
function preparedf!(
    infile::String, delim::String, datatype::String, idcol::String, firstcolpos::Int,
    testfraction::Float64=1.0, randseed=1892
)
    @info "$(loggingtime())\tpreparedf" infile delim datatype idcol firstcolpos testfraction randseed

    # read the file into a df
    if datatype == "Reads"
        df = DataFrame(CSV.File(infile, delim=delim))
    elseif datatype == "Proteins"
        # note: "Protein" is actually the ID of a *unique* protein!!!
        df1 = DataFrame(CSV.File(infile, delim=delim, select=collect(1:firstcolpos-1), types=Dict("Protein" => String, "Reads" => String)))  
        df1[!, "Protein"] = InlineString.(df1[!, :Protein])
        # make sure columns of AAs containing only Ts aren't parsed as boolean columns
        df2 = DataFrame(CSV.File(infile, delim=delim, drop=collect(1:firstcolpos-1), types=String))
        df = hcat(df1, df2)
    else
        throw("`datatype` must be either `Reads` or `Proteins`")
    end

    # take a subset of the df for testing purposes
    if testfraction < 1.0
        # Random.seed!(randseed)
        nrows = size(df, 1)
        nsamplerows = Int(round(testfraction * nrows))
        # samplerows = sample(1:nrows, nsamplerows, replace=false)
        samplerows = sample(MersenneTwister(randseed), 1:nrows, nsamplerows, replace=false)
        df = df[samplerows, :]
    end


    # flatten the df by exploding the `Reads` col, 
    # which denotes the reads supporting the unique observation (read / unique) the row represents
    transform!(df, :Reads => (x -> split.(x, ",")) => :Reads)
    df = flatten(df, :Reads)
    if datatype == "Reads"
        df = hcat(select(df, idcol), df[:, firstcolpos:end])
        firstcolpos = 2
    elseif datatype == "Proteins"
        
        # df = hcat(select(df, idcol), toAAset.(df[:, firstcolpos:end]))
        # firstcolpos = 2
        
        # include also the name of every read supporting a unique protein
        rename!(df, "Reads" => "Read")
        df = hcat(select(df, idcol, "Read"), toAAset.(df[:, firstcolpos:end]))
        firstcolpos = 3

    else
        throw("datatype must be either `Reads` or `Proteins`")
    end

    # remove uniformative cols 
    # (altough we also apply the function on the idcol it shouldn't matter for the idcol, 
    # if it has more than one unique value)
    informativedf = df[:, map(col -> length(unique(col)) > 1, eachcol(df))]
    if size(informativedf)[1] == 0
        df = df[1, :]
    end

    return df, firstcolpos
end