# Linopy2Plasmo.jl

Converts [linopy](https://linopy.readthedocs.io/) optimization models stored as NetCDF files into [Plasmo.jl](https://github.com/plasmo-dev/Plasmo.jl) graph-based optimization problems. This enables the use of graph-decomposition algorithms such as Benders decomposition on models originally built with linopy.

## Installation

### 1. Julia dependencies

Download the package and import it:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/leonardgoeke/Linopy2Plasmo.jl.git"))
Pkg.instantiate()
```

## Usage

Run the code described below from the Julia REPL after activating the project environment.

### 1. Load the packages

Load the main packages installed above:

```julia
using Linopy2Plasmo
```

Load other package used for the testing here:

```julia
using Plasmo, PlasmoBenders, PlasmoPlots
using Gurobi, CSV, DataFrames
```

### 2. Load a linopy model from NetCDF

As a test, `lin2plasObj` reads the NetCDF file provided in the data folder, extracts all variables, constraints, and the objective function, and maps set indices to integer keys for efficient graph construction.

```julia
lin2plas_obj = lin2plasObj("data/testProblem.nc")
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
structureIntoNodes!(lin2plas_obj, split_tup)
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
def_dic[:top] = filter(x -> x[1] in (0,) || x[2] == 0, lin2plas_obj.nodes) # all nodes with 0 for for the first or second set form the top-problem
for (i, j) in enumerate([1,2,3]) # creates three sub-problems depending on value of first node
    def_dic[Symbol(:sub,i)] = filter(x -> x[1] == j && x[2] != 0, lin2plas_obj.nodes)
end

```

### 5. Build the Plasmo optimization problem

`createOptProblem!` creates the `OptiGraph`, populates it with `OptiNode`s, adds variables, constraints, linking constraints, and node objectives.

```julia
createOptProblem!(lin2plas_obj, def_dic)
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
var_df[!,:value] = map(x -> Plasmo.value(benders_opt, lin2plas_obj.varMap[(x.key, x.subGraph[1])]), eachrow(var_df))
CSV.write("capacity_results.csv", select!(var_df, Not(:subGraph)))
```

## API reference

### Types

| Type | Description |
|---|---|
| `lin2plasObj` | Read-in of a linopy NetCDF model and to be converted into a Plasmo graph |
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
- `structureIntoNodes!` is the main performance bottleneck and should be optimized for large models.
