# Linopy2Plasmo.jl

Converts [linopy](https://linopy.readthedocs.io/) optimization models stored as NetCDF files into [Plasmo.jl](https://github.com/lanl-ansi/Plasmo.jl) graph-based optimization problems. This enables the use of graph-decomposition algorithms such as Benders decomposition on models originally built with linopy.

## Package structure

```
linopy2plasmo/
├── src/
│   ├── Linopy2Plasmo.jl    # module entry point
│   ├── objects.jl          # cnsObj and lin2plasObj data structures
│   └── functions.jl        # conversion and graph-building functions
├── scripts/
│   └── convert.jl          # standalone conversion script (CLI)
├── Project.toml
├── Manifest.toml
└── CondaPkg.toml           # Python environment (linopy)
```

> **Note:** The example NetCDF file (`PI_small.nc`) is not included in this repository due to its size (~2.4 GB). Provide your own linopy model exported with `model.to_netcdf("model.nc")`.

## Installation

### 1. Julia dependencies

Activate the project environment and instantiate it:

```julia
using Pkg
Pkg.activate("path/to/linopy2plasmo")
Pkg.instantiate()
```

### 2. Python dependency (linopy)

The package uses [CondaPkg.jl](https://github.com/cjdoris/CondaPkg.jl) to manage the Python environment. After activating the project, run:

```julia
using CondaPkg
CondaPkg.resolve()
```

## Usage

The workflow below follows [`scripts/convert.jl`](scripts/convert.jl). Run it from the Julia REPL after activating the project environment.

### 1. Load the package

```julia
using Linopy2Plasmo
using Plasmo, PlasmoBenders, PlasmoPlots
using Gurobi, CSV
```

### 2. Load a linopy model from NetCDF

`lin2plasObj` reads the NetCDF file, extracts all variables, constraints, and the objective function, and maps set indices to integer keys for efficient graph construction.

```julia
lin2plas_obj = lin2plasObj("PI_small.nc")
```

You can inspect the extracted data directly:

```julia
lin2plas_obj.sets[:set_technologies]          # set name → integer index
lin2plas_obj.var[:capacity]                   # variable data as DataFrame
lin2plas_obj.cns[:constraint_couple_storage_level]  # constraint data
```

### 3. Structure the model into graph nodes

`structureIntoNodes!` assigns every variable and constraint to one or more nodes in the Plasmo graph, based on a tuple of set names that defines the node hierarchy. Each unique combination of set indices becomes a node.

```julia
split_tup = (:set_time_steps_yearly, :set_nodes, :set_time_steps_operation, :set_time_steps_storage)
@time structureIntoNodes!(lin2plas_obj, split_tup)
```

After this call, `lin2plas_obj.nodes` contains all node tuples and `lin2plas_obj.varNode` maps each variable ID to the nodes it appears in.

### 4. Define the subgraph layout

Subgraphs group nodes for decomposition. The simplest option puts everything in one subgraph:

```julia
def_dic = Dict(:all => lin2plas_obj.nodes)
```

For Benders decomposition you need at least two subgraphs — a master (`:top`) and one or more sub-problems. The example below creates one sub-problem per year × location combination:

```julia
def_dic = Dict{Symbol, Vector{Tuple{Vararg{Int}}}}()
def_dic[:top] = filter(x -> x[2] == 0 || (x[3] == 0 && x[4] == 0), lin2plas_obj.nodes)
nonTop_arr    = filter(x -> !(x[2] == 0 || (x[3] == 0 && x[4] == 0)), lin2plas_obj.nodes)

for (idx, (i, j)) in enumerate(unique(getindex.(nonTop_arr, 1)))
    def_dic[Symbol(:sub, idx)] = filter(x -> x[1] == i && x[2] == j, nonTop_arr)
end
```

### 5. Build the Plasmo optimization problem

`createOptProblem!` creates the `OptiGraph`, populates it with `OptiNode`s, adds variables, constraints, linking constraints, and node objectives.

```julia
@time createOptProblem!(lin2plas_obj, def_dic)
```

Optional — visualize the graph structure:

```julia
matrix_plot(lin2plas_obj.mainGraph, subgraph_colors = true)
```

### 6. Solve

**Standard solve** (single solver across all subgraphs):

```julia
foreach(x -> set_optimizer(lin2plas_obj.subGraphs[x], Gurobi.Optimizer), keys(lin2plas_obj.subGraphs))
set_optimizer(lin2plas_obj.mainGraph, Gurobi.Optimizer)
optimize!(lin2plas_obj.mainGraph)
```

**Benders decomposition** (requires a `:top` / `:sub` subgraph layout):

```julia
benders_opt = BendersAlgorithm(
    lin2plas_obj.mainGraph, lin2plas_obj.subGraphs[:top];
    solver     = optimizer_with_attributes(Gurobi.Optimizer),
    add_slacks = true,
    max_iters  = 1000,
    regularize = true,
)
run_algorithm!(benders_opt)
```

### 7. Extract and write results

`replaceSetColumns` converts the integer set indices stored internally back to their original string labels.

```julia
var_df = replaceSetColumns(copy(lin2plas_obj.var[:capacity]), lin2plas_obj.revSets)

# standard solve
var_df[!, :value] = map(
    x -> Plasmo.value(lin2plas_obj.varMap[(x.key, x.subGraph[1])]),
    eachrow(var_df),
)

# Benders solve (pass the algorithm object instead)
# var_df[!, :value] = map(
#     x -> Plasmo.value(benders_opt, lin2plas_obj.varMap[(x.key, x.subGraph[1])]),
#     eachrow(var_df),
# )

CSV.write("capacity_results.csv", var_df)
```

## API reference

### Types

| Type | Description |
|---|---|
| `lin2plasObj(path)` | Reads a linopy NetCDF model and returns the conversion object |
| `cnsObj` | Internal representation of a single constraint (variables, RHS, sense) |

### Key fields of `lin2plasObj`

| Field | Type | Description |
|---|---|---|
| `sets` | `Dict{Symbol, Dict{String,Int}}` | Set name → element → integer index |
| `revSets` | `Dict{Symbol, Dict{Int,String}}` | Reverse mapping (index → element name) |
| `var` | `Dict{Symbol, DataFrame}` | Variable data per variable name |
| `cns` | `Dict{Symbol, DataFrame}` | Constraint data per constraint name |
| `nodes` | `Array{Tuple}` | All node tuples (one per unique set combination) |
| `varNode` | `Dict{Int, Array{Tuple}}` | Variable ID → list of nodes it appears in |
| `varMap` | `Dict{(Int,Tuple), VariableRef}` | (variable ID, node) → JuMP variable reference |
| `mainGraph` | `OptiGraph` | The top-level Plasmo graph |
| `subGraphs` | `Dict{Symbol, OptiGraph}` | Named subgraphs |
| `linkCns` | `DataFrame` | All linking constraints with their node pairs |

### Functions

| Function | Description |
|---|---|
| `structureIntoNodes!(obj, split_tup)` | Assigns variables and constraints to nodes based on the set hierarchy |
| `createOptProblem!(obj, def_dic; wrtStr, isPlasmo)` | Builds the full Plasmo OptiGraph |
| `replaceSetColumns(df, revSets)` | Converts integer indices in a DataFrame back to set-element strings |
| `sortNodes(nodes)` | Groups node tuples by hierarchy level (number of zero entries) |
| `checkISS(obj)` | Prints IIS constraints for infeasible JuMP (non-Plasmo) models |

## Known limitations

- Binary and integer variables are not yet supported (continuous variables only).
- `structureIntoNodes!` is the main performance bottleneck (~95 s / 1.6 GB allocations on an 11 k-constraint model).
