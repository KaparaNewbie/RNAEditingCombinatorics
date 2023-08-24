#= 

G = (V, E)
V = {a, b, c, d}
E = {(b, d), (c, d)}

Acceptable solutions:
1. {a, b}
2. {a, c}
3. {b, c}
4. {a, b, c} (which is also the best solution)

Packages: JuMP, Clp
=# 



# using JuMP
# # using Clp
# using HiGHS

# # model = Model(Clp.Optimizer)
# model = Model(HiGHS.Optimizer)


# @variable(model, 0 <= a <= 1, Int)
# @variable(model, 0 <= b <= 1, Int)
# @variable(model, 0 <= c <= 1, Int)
# @variable(model, 0 <= d <= 1, Int)

# @objective(model, Max, a + b + c + d)

# @constraint(model, c1, 0 <= b + d <= 1)
# @constraint(model, c2, 0 <= c + d <= 1)

# print(model)

# optimize!(model)

# termination_status(model)

# objective_value(model)

# value(a)
# value(b)
# value(c)
# value(d)



using JuMP
using HiGHS


model = Model(HiGHS.Optimizer)


sort_edge_by_vertices_names(u, v) = u < v ? (u, v) : (v, u)

G = Dict(
    "a" => Set(),
    "b" => Set(["d"]),
    "c" => Set(["d"]),
    "d" => Set(["b", "c"])
)
V = keys(G)
E = Set([sort_edge_by_vertices_names(u, v) for u ∈ V for v ∈ G[u] if u != v])

@variable(model, x[V], Bin)

@objective(model, Max, sum(x))

for (u, v) ∈ E
    @constraint(model, 0 <= x[u] + x[v] <= 1)
end

print(model)

optimize!(model)

termination_status(model)

objective_value(model)

for u in V
    println(u, " = ", value(x[u]))
end
