# using Transducers  # for tcollect
using Statistics
using Random
using CairoMakie
CairoMakie.activate!()

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

# meaneditingfreq(M) = mean(sum(M, dims=2) / size(M, 2))
# calctotaleditingfreq(M) == meaneditingfreq(M)

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
            # println("i = $i, j = $j, d = $d")
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
    # editingfreq = meaneditingfreq(M)
    # step = 0
    while length(remainingreads) > 0 && totaleditfreq < maxeditfreq
    # while length(remainingreads) > 0 && editingfreq < maxeditfreq
        # step = step + 1
        Threads.@threads for i ∈ remainingreads
            row = view(M, i, :)
            editrow!(row, editingprob)
        end
        remainingreads = removereads(remainingreads, readremovalprob)
        totaleditfreq = calctotaleditingfreq(M)
        # editingfreq = meaneditingfreq(M)
    end

    cors = cor(M)
    return cors
end

function simulate(N, K, editingprob, readremovalprob, mask = true)
    M = simulatereads(N, K)
    remainingreads = collect(1:N)
    # totaleditfreq = calctotaleditingfreq(M)
    # editingfreq = meaneditingfreq(M)
    # step = 0
    while length(remainingreads) > 0
    # while length(remainingreads) > 0 && editingfreq < maxeditfreq
        # step = step + 1
        Threads.@threads for i ∈ remainingreads
            row = view(M, i, :)
            editrow!(row, editingprob)
        end
        remainingreads = removereads(remainingreads, readremovalprob)
        # totaleditfreq = calctotaleditingfreq(M)
        # editingfreq = meaneditingfreq(M)
    end

    cors = cor(M)
    # if mask
    #     cors = mask_tril_with_nan(cors)
    # end
    return cors
end



settitle(editingprob, readremovalprob) = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%"
settitle(editingprob, readremovalprob, i, j) = "$(settitle(editingprob, readremovalprob))\ni = $i, j = $j"
settitle(editingprob, readremovalprob, expectedcor::Float64, ndigits::Int=3) = "$(settitle(editingprob, readremovalprob))\nexpected cor = $(round(expectedcor; digits=ndigits))"


intifpossible(x) = isinteger(x) ? Int(x) : x

# editingprob = 0.001
# editingprob = 100 * editingprob
# if isinteger(editingprob)
#     editingprob = Int(editingprob)
# end
# expected = expectedcor(readremovalprob, editingprob)
# round100*editingprob
# title = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%\nexpected cor = $(round(expectedcor; digits=ndigits))"
# title1 = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%"

function settitle3(editingprob, readremovalprob)
    editingprob = intifpossible(100 * editingprob)
    readremovalprob = intifpossible(100 * readremovalprob)
    title = "edit site = $(editingprob)%\nremove read = $(readremovalprob)%"
end


function settitle3(editingprob, readremovalprob, expectedcor::Float64, ndigits::Int=3)
    title = "$(settitle2(editingprob, readremovalprob))\nexpected cor = $(round(expectedcor; digits=ndigits))"
end



# expectedcor(readremovalprob, editingprob) = editingprob / (readremovalprob + 2 * editingprob)
expectedcor(readremovalprob, editingprob) =  editingprob / (2 * editingprob + readremovalprob)
expectedcor(0.003, 0.001)
# 3/7


N = 10_000 # number of reads
K = 10 # number of editing sites per read
# N, K = 10, 3



editingprobs = [0.01, 0.02, 0.03] # chance of a site being edited at each step
readremovalprobs = [0.01, 0.02, 0.03] # chance of a read being removed from the simulation at each step
# editingprobs = [0.001, 0.002, 0.003] # chance of a site being edited at each step
# readremovalprobs = [0.001, 0.002, 0.003] # chance of a read being removed from the simulation at each step
# maxeditfreqs = [0.1, 0.5, 0.9] # maximum fraction of edited sites across all reads and sites


# corsdict = Dict(
#     (editingprob, readremovalprob, maxeditfreq) => simulate(N, K, editingprob, readremovalprob, maxeditfreq) 
#     for editingprob in editingprobs, readremovalprob in readremovalprobs, maxeditfreq in maxeditfreqs
# )

corsdict = Dict(
    (editingprob, readremovalprob) => simulate(N, K, editingprob, readremovalprob) 
    for editingprob in editingprobs, readremovalprob in readremovalprobs
)
# corsdict[(0.01, 0.03)]
# corsdict[(0.03, 0.01)]
maskedcorsdict = Dict(
    (editingprob, readremovalprob) => mask_tril_with_nan(corsdict[(editingprob, readremovalprob)])
    for (editingprob, readremovalprob) in keys(corsdict)
)


abs_min = minimum(x -> minimum(filter(!isnan, x)), values(maskedcorsdict))
abs_max = maximum(x -> maximum(filter(!isnan, x)), values(maskedcorsdict))
# abs_min = floor(abs_min, digits = 1)
# abs_max = ceil(abs_max, digits = 1)
joint_limits = (abs_min, abs_max)

xs = ys = collect(1:K)

# colormap = cgrad(:Spectral, 20, categorical = true)
colormap = cgrad(:Spectral,)

fig = Figure(size = (600, 600))
for (i, editingprob) in enumerate(editingprobs)
    xaxisdetails = i == length(editingprobs)
    for (j, readremovalprob) in enumerate(readremovalprobs)
        z = maskedcorsdict[(editingprob, readremovalprob)]
        title = settitle(editingprob, readremovalprob)
        yaxisdetails = j == 1
        ax = Axis(
            fig[i, j], 
            subtitle = title, 
            xticksvisible = xaxisdetails, xticklabelsvisible = xaxisdetails,
            yticksvisible = yaxisdetails, yticklabelsvisible = yaxisdetails
        )
        heatmap!(ax, xs, ys, z,  colorrange = joint_limits, colormap = colormap)
        hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
    end
end
Label(fig[end+1, :], "Site")  # a axis title
Label(fig[begin:end-1, 0], "Site", rotation = pi/2)  # y axis title
Colorbar(
    fig[begin:end-1, end+1], 
    limits = joint_limits, 
    label = "Pearson's r", 
    colormap = colormap
)
fig






editingprobs = [0.001, 0.002, 0.003] # chance of a site being edited at each step
readremovalprobs = [0.001, 0.002, 0.003] # chance of a read being removed from the simulation at each step
corsdict = Dict(
    (editingprob, readremovalprob) => simulate(N, K, editingprob, readremovalprob) 
    for editingprob in editingprobs, readremovalprob in readremovalprobs
)
maskedcorsdict = Dict(
    (editingprob, readremovalprob) => mask_tril_with_nan(corsdict[(editingprob, readremovalprob)])
    for (editingprob, readremovalprob) in keys(corsdict)
)




fig = Figure(size = (600, 600))
abs_min = minimum(x -> minimum(filter(!isnan, x)), values(maskedcorsdict))
abs_max = maximum(x -> maximum(filter(!isnan, x)), values(maskedcorsdict))
# abs_min_floor = floor(abs_min, digits = 1)
# abs_max_ceil = ceil(abs_max, digits = 1)
# limits = (abs_min_floor, abs_max_ceil, nothing, nothing)
limits = (abs_min, abs_max, nothing, nothing)
axes = []
for (i, editingprob) in enumerate(editingprobs)
    # xaxisdetails = i == length(editingprobs)
    # rowslabel = "edit site = $(intifpossible(100 * editingprob))%"
    for (j, readremovalprob) in enumerate(readremovalprobs)
        yaxisdetails = j == 1
        z = maskedcorsdict[(editingprob, readremovalprob)]
        z = filter(!isnan, z)
        # title = settitle(editingprob, readremovalprob)
        expected = expectedcor(readremovalprob, editingprob)
        title = settitle3(editingprob, readremovalprob, expected, 2)
        # title = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%\nexpected cor = $(round(expectedcor; digits=ndigits))"
        # title = "expected cor = $(round(expected; digits=2))"
        # yaxisdetails = j == 1
        ax = Axis(
            fig[i, j], 
            subtitle = title, 
            limits = limits,
            # xticksvisible = xaxisdetails, xticklabelsvisible = xaxisdetails,
            # yticksvisible = yaxisdetails, yticklabelsvisible = yaxisdetails
        )
        push!(axes, ax)
        hist!(
            ax, z; 
            # label = title
        )
        # axislegend(ax; position = :rt)
        # hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
        vlines([1, 2, 3])
        
    end
end
linkxaxes!(axes...)
linkyaxes!(axes...)
Label(fig[end+1, :], "Pearson's r")  # a axis title
Label(fig[begin:end-1, 0], "Sites", rotation = pi/2)  # y axis title
fig


# As = Bs = [1, 2]
# Xs = Ys = [1, 2]

# cors_dict = Dict(
#     (1, 1) => [0.3 0.3; 0.3 0.3],
#     (1, 2) => [0.2 0.2; 0.2 0.2],
#     (2, 1) => [0.6 0.6; 0.6 0.6],
#     (2, 2) => [0.3 0.3; 0.3 0.3]
# )

# joint_limits = (0.1, 0.6)

# colormap = cgrad(:Spectral, 20, categorical = true)

# fig = Figure(size = (400, 400))

# for a in As
#     for b in Bs 
#         zs = cors_dict[(a, b)]
#         println("a = $a, b = $b")
#         println(zs)
#         title = "a = $a\nb = $b"
#         ax = Axis(fig[a, b], title = title)
#         # heatmap!(ax, Xs, Ys, zs,  colorrange = joint_limits)
#         heatmap!(ax, Xs, Ys, zs,  colorrange = joint_limits, colormap = colormap)

#     end
# end

# Colorbar(
#     fig[:, 3], 
#     limits = joint_limits, 
#     label = "Pearson's r", 
#     # colormap = cgrad(:Spectral, 20, categorical = true)
#     colormap = colormap
# )

# fig









# Label(fig[0, :], "maxeditfreq = $maxeditfreq", fontsize = 20)




# figs = []
# for maxeditfreq in maxeditfreqs
#     fig = Figure(size = (600, 600))
#     i = j = 1
#     editingprob = editingprobs[i]
#     readremovalprob = readremovalprobs[j]
#     z = corsdict[(editingprob, readremovalprob, maxeditfreq)]
#     z = mask_tril_with_nan(z)
#     title = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%"
#     ax = Axis(fig[i, j], title = title)
#     heatmap!(ax, xs, ys, z,  colorrange = joint_limits)
#     hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
#     for (i, editingprob) in enumerate(editingprobs)
#         for (j, readremovalprob) in enumerate(readremovalprobs)
#             i == 1 && j == 1 && continue
#             z = corsdict[(editingprob, readremovalprob, maxeditfreq)]
#             z = mask_tril_with_nan(z)
#             title = "edit site = $(Int(100*editingprob))%\nremove read = $(Int(100*readremovalprob))%"
#             ax = Axis(fig[i, j], title = title)
#             heatmap!(ax, xs, ys, z,  colorrange = joint_limits)
#             hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
#         end
#     end
#     Colorbar(fig[:, end+1], limits = joint_limits, label = "Pearson's r", colormap = cgrad(:Spectral, 20, categorical = true))
#     Label(fig[0, :], "maxeditfreq = $maxeditfreq", fontsize = 20)
#     push!(figs, fig)
# end

# figs[1]
# figs[2]
# figs[3]



