# basic slack versiob

# pmaps (didn't actually used these llinks in here)
# https://discourse.julialang.org/t/pmap-in-function-module/25896
# https://stackoverflow.com/questions/44741667/julia-equivalent-of-python-multiprocessing-pool-map


using Distributed
addprocs(4)

using Transducers
using Random
using StatsBase
# using BenchmarkTools


Random.seed!(1892) # for reproducibility
M = sample([-1, 0, 1], (1000, 10)) # a matrix with 1000 rows, 10 columns
k = "some_key"
d = Dict(k => 3)

sumrow(row, d, k) = sum(row) * d[k] # an arbitrary function that accepts a number of arguments


s1 = dcollect(sumrow(row, d, k) for row in eachrow(M))


s2 = dcollect(sumrow(row, d, k) for row in eachrow(M))






# JULIA_PROJECT=. julia

# JULIA_PROJECT=. julia -t 20



"""
using Distributed
addprocs(20)  # or addprocs(8; exeflags="--project")
@everywhere using Transducers
using Random

nrows, mcols = 10_000, 10_000
Random.seed!(1892) # for reproducibility
M = rand(Int, (nrows, mcols))

@everywhere begin
    k = "some_key"
    d = Dict(k => 100)
    sumrow(row, d, k) = sum((row .^ 3) .+ (d[k]^2 .* row .+ row)) # an arbitrary function that accepts a number of arguments
    sumrow2(row, x) = sum((row .^ 3) .+ (x ^ 2 .* row .+ row)) # an arbitrary function that accepts a number of arguments
end

s1 = [sumrow(row, d, k) for row in eachrow(M)] # sequential
s2 = dcollect(sumrow(row, d, k) for row in eachrow(M)) # distributed
s3 = tcollect(sumrow(row, d, k) for row in eachrow(M)) # threaded
@assert s1 == s2 == s3

using BenchmarkTools
@benchmark [sumrow(row, d, k) for row in eachrow(M)]  # sequential
@benchmark dcollect(sumrow(row, d, k) for row in eachrow(M))  # distributed
@benchmark tcollect(sumrow(row, d, k) for row in eachrow(M))  # threaded

@benchmark [sumrow2(row, d[k]) for row in eachrow(M)]  # sequential
@benchmark dcollect(sumrow2(row, d[k]) for row in eachrow(M))  # distributed
@benchmark tcollect(sumrow2(row, d[k]) for row in eachrow(M))  # threaded
"""


# everywhere version



using Distributed
addprocs(4)

@everywhere begin
    using Transducers
    using Random
    using StatsBase
    # using BenchmarkTools


    Random.seed!(1892) # for reproducibility
    M = sample([-1, 0, 1], (1000, 10)) # a matrix with 1000 rows, 10 columns
    k = "some_key"
    d = Dict(k => 3)

    sumrow(row, d, k) = sum(row .* d[k]^2 .* row)  # an arbitrary function that accepts a number of arguments
end

s1 = dcollect(sumrow(row, d, k) for row in eachrow(M))


s2 = dcollect(sumrow(row, d, k) for row in eachrow(M))



# a = 43
# Random.seed!(1892) # for reproducibility
# M = sample([-1, 0, 1], (1000, 10)) # a matrix with 1000 rows, 10 columns

# sumrow(row, firstcolpos) = sum(row[firstcolpos:end]) # an arbitrary function that accepts a number of arguments

# addprocs(8)
# firstcolpos = 2
# s = dcollect(sumrow(row, firstcolpos) for row in eachrow(M))





# s2 = [sumrow(row, firstcolpos) for row in eachrow(M)]



# s == s2

# @benchmark dcollect(sumrow(row, firstcolpos) for row in eachrow(M))
# @benchmark [sumrow(row, firstcolpos) for row in eachrow(M)]



# @everywhere sumrow2(row) = sum(row)
# pmap(x -> sumrow(x), (x for x in c))
# s2 = pmap(x -> sumrow(x, 1), c)
# s2 = pmap(x -> sumrow(x, 1), (x for x in c))





