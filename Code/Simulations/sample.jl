

function get_graph_sample(G::Dict, fraction::Float64, nsamplerows::Int64, df::DataFrame, idcol::String)
	@info "$(loggingtime())\tget_graph_sample" fraction myid()
	if fraction < 1.0
		sampleG = get_graph_sample(G, nsamplerows, df, idcol)
	else  # fraction == 1.0, so there's no need to create a sub-graph
		sampleG = G
	end
	return sampleG
end


function get_graph_sample(G::Dict, nsamplerows::Int64, df::DataFrame, idcol::String)
	# sample a fraction of the rows (non-unique reads/proteins)
	samplerows = sample(collect(1:size(df, 1)), nsamplerows, replace = false)
	# get their corresponding unique ids
	sampleids = ThreadsX.unique(df[samplerows, idcol])
	# also get all corresponding reads' names
	samplereads = df[samplerows, "Read"]
	# assemble sub neighborhood lists of uncompatible unique sampled rows by using the pre-computed complete graph
	sampleG = subgraph(G, sampleids)
	return sampleG
end


function get_graph_sample_and_available_reads(G::Dict, fraction::Float64, nsamplerows::Int64, df::DataFrame, idcol::String)
	@info "$(loggingtime())\tget_graph_sample_and_available_reads" fraction myid()
	if fraction < 1.0
		sampleG, availablereads = get_graph_sample_and_available_reads(G, nsamplerows, df, idcol)
	else  # fraction == 1.0, so there's no need to create a sub-graph
		sampleG = G
		# print(df)
		availablereads = df[:, "Read"]
	end
	return sampleG, availablereads
end


function get_graph_sample_and_available_reads(G::Dict, nsamplerows::Int64, df::DataFrame, idcol::String)
	# sample a fraction of the rows (non-unique reads/proteins)
	samplerows = sample(collect(1:size(df, 1)), nsamplerows, replace = false)
	# get their corresponding unique ids
	sampleids = ThreadsX.unique(df[samplerows, idcol])
	# also get all corresponding avaiable reads' names
	# print(df)
	availablereads = df[samplerows, "Read"]
	# assemble sub neighborhood lists of uncompatible unique sampled rows by using the pre-computed complete graph
	sampleG = subgraph(G, sampleids)
	return sampleG, availablereads
end







function get_graph_sample_and_available_reads(
	G::Dict, fraction::Float64, samplerows::Vector{Int64}, df::DataFrame, idcol::String,
)
	@info "$(loggingtime())\tget_graph_sample_and_available_reads" fraction myid()
	if fraction < 1.0
		sampleG, availablereads = get_graph_sample_and_available_reads(G, samplerows, df, idcol)
	else  # fraction == 1.0, so there's no need to create a sub-graph
		sampleG = G
		# print(df)
		availablereads = df[:, "Read"]
	end
	return sampleG, availablereads
end

function get_graph_sample_and_available_reads(
	G::Dict, samplerows::Vector{Int64}, df::DataFrame, idcol::String,
)
	# # make sure the fraction of the required rows (non-unique reads/proteins) 
	# # is consistent with the size of the input
	# @assert length(samplerows) ≈ fraction * size(df, 1)

	# get the corresponding unique ids of input samplerows
	sampleids = ThreadsX.unique(df[samplerows, idcol])
	# also get all corresponding avaiable reads' names
	# print(df)
	availablereads = df[samplerows, "Read"]
	# assemble sub neighborhood lists of uncompatible unique sampled rows by using the pre-computed complete graph
	sampleG = subgraph(G, sampleids)
	return sampleG, availablereads
end




"""
	subgraph(G, V2)

Return a subgraph `G2 ⊆ G` induced by the vertices `V2 ⊆ V`, where `V` are the original vertices of `G`.
"""
function subgraph(G, V2)
	V2eltype = eltype(V2)
	# vertices in G2
	G2 = Dict(u => Set{V2eltype}() for u ∈ V2)
	# edges in G2
	length(G2) == length(V2) || error("The vertices in V2 should be unique.")
	Threads.@threads for u ∈ V2
		Threads.@threads for v ∈ V2 # compare vertices u & v if u != v
			u != v && v ∈ G[u] && push!(G2[u], v)
		end
	end
	return G2
end
