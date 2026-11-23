# Exact scoped Figure S19 producer. Run with Julia 1.7.2 and --threads=1.
# Scientific functions and comprehension order come from
# unifinished_editing_isoforms.jl; plotting and unrelated experiments are omitted.
using Random, Statistics, Serialization, LinearAlgebra, StatsBase
@assert VERSION == v"1.7.2"
@assert Threads.nthreads() == 1
@assert length(ARGS) == 2 "Supply new table and intermediate directories"
const OUTPUT = abspath(ARGS[1])
@assert !ispath(OUTPUT) "Output directory must be new; approved results are never overwritten"
mkpath(OUTPUT)
const WORK_OUTPUT = abspath(ARGS[2])
@assert !ispath(WORK_OUTPUT) "Intermediate directory must be new"
mkpath(WORK_OUTPUT)

function simulate_and_save(N, K, editingprob, readremovalprob)
    key = "$(editingprob)_$(readremovalprob)"
    serialize(joinpath(WORK_OUTPUT, "rng_before_$(key).jls"), copy(Random.default_rng()))
    started = time()
    matrix = simulate(N, K, editingprob, readremovalprob)
    serialize(joinpath(WORK_OUTPUT, "correlations_$(key).jls"), matrix)
    serialize(joinpath(WORK_OUTPUT, "rng_after_$(key).jls"), copy(Random.default_rng()))
    println("Completed ",key," in ",round(time()-started; digits=2)," seconds")
    flush(stdout)
    return matrix
end

Random.seed!(1892)

"""Create a matrix of `N` reads with `K` editing sites, all zeroed out (unedited) at first."""
simulatereads(N, K) = zeros(Bool, N, K)

readremains(readremovalprob) = rand() > readremovalprob

"""
Remove reads from the vector `remainingreads` with probability `readremovalprob`.
"""
function removereads(remainingreads, readremovalprob)
	remainingindices = [readremains(readremovalprob) for _ ∈ eachindex(remainingreads)]
	return remainingreads[remainingindices]
end

function calctotaleditingfreq(M)
	return sum(M) / length(M)
end

edit(editingprob) = rand() <= editingprob

function editrow!(row, editingprob)
	editedrow = [edit(editingprob) for _ ∈ eachindex(row)]
	setindex!(row, row .| editedrow, eachindex(row))
end

"""
Mask the lower triangular part of a matrix `A` with NaNs.
"""
function mask_tril_with_nan(A)
	Avals = []
	d = 1
	for i in 1:size(A, 1)
		for j in 1:size(A, 2)
			x = A[i, j]
			if j >= i
				x = NaN
			end
			push!(Avals, x)
		end
		d += 1
	end
	return reshape(Avals, size(A, 1), size(A, 2))
end

function simulate(N, K, editingprob, readremovalprob, maxeditfreq)
	M = simulatereads(N, K)
	remainingreads = collect(1:N)
	totaleditfreq = calctotaleditingfreq(M)
	while length(remainingreads) > 0 && totaleditfreq < maxeditfreq
		Threads.@threads for i ∈ remainingreads
			row = view(M, i, :)
			editrow!(row, editingprob)
		end
		remainingreads = removereads(remainingreads, readremovalprob)
		totaleditfreq = calctotaleditingfreq(M)
	end

	serialize(joinpath(WORK_OUTPUT, "reads_$(editingprob)_$(readremovalprob).jls"), M)
	cors = cor(M)
	return cors
end

function simulate(N, K, editingprob, readremovalprob, mask = true)
	M = simulatereads(N, K)
	remainingreads = collect(1:N)
	while length(remainingreads) > 0
		Threads.@threads for i ∈ remainingreads
			row = view(M, i, :)
			editrow!(row, editingprob)
		end
		remainingreads = removereads(remainingreads, readremovalprob)
	end

	serialize(joinpath(WORK_OUTPUT, "reads_$(editingprob)_$(readremovalprob).jls"), M)
	cors = cor(M)
	return cors
end


expectedcor(readremovalprob, editingprob) = editingprob / (2 * editingprob + readremovalprob)
N = 10_000
K = 10
editingprobs = [0.001, 0.002, 0.003] # chance of a site being edited at each step
readremovalprobs = [0.001, 0.002, 0.003] # chance of a read being removed from the simulation at each step
corsdict = Dict(
	(editingprob, readremovalprob) => simulate_and_save(N, K, editingprob, readremovalprob)
	for editingprob in editingprobs, readremovalprob in readremovalprobs
)
maskedcorsdict = Dict(
	(editingprob, readremovalprob) => mask_tril_with_nan(corsdict[(editingprob, readremovalprob)])
	for (editingprob, readremovalprob) in keys(corsdict)
)



function pick_hist_edges(vals, bins)
    if bins isa Int
        mi, ma = float.(extrema(vals))
        if mi == ma
            return (mi - 0.5):(ma + 0.5)
        end
        # hist is right-open, so to include the upper data point, make the last bin a tiny bit bigger
        ma = nextfloat(ma)
        return range(mi, ma, length = bins+1)
    else
        if !issorted(bins)
            error("Histogram bins are not sorted: $bins")
        end
        return bins
    end
end

open(joinpath(OUTPUT,"s19_correlations.tsv"),"w") do io
    println(io,"EditingProbability\tReadRemovalProbability\tRow\tColumn\tCorrelation\tDisplayed\tExpectedCorrelation")
    for editingprob in editingprobs, removalprob in readremovalprobs
        matrix=corsdict[(editingprob,removalprob)]
        masked=maskedcorsdict[(editingprob,removalprob)]
        @assert maximum(abs.(matrix-transpose(matrix))) < 1e-12
        @assert all(isfinite,matrix) && all(diag(matrix).≈1)
        @assert count(!isnan,masked)==45
        for col in 1:K, row in 1:K
            println(io,join((editingprob,removalprob,row,col,matrix[row,col],!isnan(masked[row,col]),expectedcor(removalprob,editingprob)),'\t'))
        end
    end
end
open(joinpath(OUTPUT,"s19_histograms.tsv"),"w") do io
    println(io,"EditingProbability\tReadRemovalProbability\tBin\tLeft\tRight\tCount\tExpectedCorrelation")
    for editingprob in editingprobs, removalprob in readremovalprobs
        values=Float64.(filter(!isnan,maskedcorsdict[(editingprob,removalprob)]))
        edges=pick_hist_edges(values,15)
        histogram=StatsBase.fit(StatsBase.Histogram,values,edges)
        @assert sum(histogram.weights)==45 && length(histogram.weights)==15
        for b in 1:15
            println(io,join((editingprob,removalprob,b,edges[b],edges[b+1],histogram.weights[b],expectedcor(removalprob,editingprob)),'\t'))
        end
    end
end
open(joinpath(OUTPUT,"s19_complete.tsv"),"w") do io
    println(io,"Setting\tValue\nSeed\t1892\nReads\t10000\nSites\t10\nThreads\t1\nJulia\t",VERSION,"\nSimulationPanels\t9\nPairwiseValuesPerPanel\t45\nBinsPerPanel\t15")
end
println("S19 simulation and histogram export complete.")
