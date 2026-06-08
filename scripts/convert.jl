# Script: convert a linopy NetCDF model to a Plasmo OptiGraph and solve it.
#
# Usage (standard solve):
#   julia --project scripts/convert.jl <path.nc> [output.csv]
#
# Usage (Benders decomposition):
#   julia --project scripts/convert.jl <path.nc> [output.csv] --benders

using Linopy2Plasmo
using Plasmo, PlasmoBenders
using Gurobi, CSV

# ---------- arguments ----------
nc_path     = length(ARGS) >= 1 ? ARGS[1] : "PI_small.nc"
output_path = length(ARGS) >= 2 ? ARGS[2] : "capacity_results.csv"
use_benders = "--benders" in ARGS

# ---------- 1. load linopy model ----------
lin2plas_obj = lin2plasObj(nc_path)

# ---------- 2. structure into nodes ----------
split_tup = (:set_time_steps_yearly, :set_nodes, :set_time_steps_operation, :set_time_steps_storage)
@time structureIntoNodes!(lin2plas_obj, split_tup)

# ---------- 3. define subgraph layout ----------
# All-in-one (single subgraph):
def_dic = Dict(:all => lin2plas_obj.nodes)

# Alternatively — one subgraph per yearly period × location:
#=
def_dic = Dict{Symbol, Vector{Tuple{Vararg{Int}}}}()
def_dic[:top] = filter(x -> x[2] == 0 || (x[3] == 0 && x[4] == 0), lin2plas_obj.nodes)
nonTop_arr   = filter(x -> !(x[2] == 0 || (x[3] == 0 && x[4] == 0)), lin2plas_obj.nodes)
for (idx, (i, j)) in enumerate(unique(getindex.(nonTop_arr, 1)))
    def_dic[Symbol(:sub, idx)] = filter(x -> x[1] == i && x[2] == j, nonTop_arr)
end
=#

# ---------- 4. build optimization problem ----------
@time createOptProblem!(lin2plas_obj, def_dic)

# ---------- 5. solve ----------
if use_benders
    benders_opt = BendersAlgorithm(
        lin2plas_obj.mainGraph,
        lin2plas_obj.subGraphs[:top];
        solver    = optimizer_with_attributes(Gurobi.Optimizer),
        add_slacks = true,
        max_iters  = 1000,
        regularize = true,
    )
    run_algorithm!(benders_opt)
else
    foreach(x -> set_optimizer(lin2plas_obj.subGraphs[x], Gurobi.Optimizer), keys(lin2plas_obj.subGraphs))
    set_optimizer(lin2plas_obj.mainGraph, Gurobi.Optimizer)
    optimize!(lin2plas_obj.mainGraph)
end

# ---------- 6. extract and write results ----------
var_df = replaceSetColumns(copy(lin2plas_obj.var[:capacity]), lin2plas_obj.revSets)

if use_benders
    var_df[!, :value] = map(
        x -> Plasmo.value(benders_opt, lin2plas_obj.varMap[(x.key, x.subGraph[1])]),
        eachrow(var_df),
    )
else
    var_df[!, :value] = map(
        x -> Plasmo.value(lin2plas_obj.varMap[(x.key, x.subGraph[1])]),
        eachrow(var_df),
    )
end

CSV.write(output_path, var_df)
println("Results written to $output_path")
