
#region # * high-level function to change plasmo object

# ! structure model into nodes based on input sets
function structureIntoNodes!(plasmo_obj::lin2plasObj, split_tup::Tuple{Vararg{Symbol}})

    plasmo_obj.setHier = split_tup

    # assign constraints and variables to nodes
    allSub_arr = Array{Tuple{Vararg{Int}}}(undef, 0)
    varSub_dic = Dict{Int, Array{Tuple{Vararg{Int}}}}()

    cns_dic = plasmo_obj.cns
    for c in keys(cns_dic)
        # scan existing sets
        inSet_arr = intersect(plasmo_obj.setHier, namesSym(cns_dic[c]))
        # add subgraph info to constraints
        cns_dic[c][!,:subGraph] .= map(y -> tuple([z in inSet_arr ? getindex(y,z) : 0 for z in plasmo_obj.setHier]...), eachrow(cns_dic[c]))
        append!(allSub_arr, unique(cns_dic[c][!,:subGraph]))
        # add subgraph info to variables
        foreach(y -> addToVarSubDic!(varSub_dic, getindex.(y.cnsObj.var,1), y.subGraph), eachrow(cns_dic[c]))
    end

    # map objective variables to nodes
    plasmo_obj.obj[!,:node] .= map(x -> sortNodes(varSub_dic[x])[1][1], plasmo_obj.obj[!,:var])


    # saver results
    plasmo_obj.setHier = split_tup
    plasmo_obj.varNode = varSub_dic
    plasmo_obj.nodes = unique(allSub_arr)

end

# ! create optimization problem from object
function createOptProblem!(plasmo_obj::lin2plasObj, def_dic::Dict{Symbol, Vector{Tuple{Vararg{Int64}}}}; wrtStr = true, isPlasmo = true)

    # ! creates graphs and nodes
    # create graphs and connect
    main_graph = isPlasmo ? OptiGraph() : Model()
    subGraph_dic = Dict(x => OptiGraph() for x in keys(def_dic))
    foreach(x -> add_subgraph!(main_graph, x), collect(values(subGraph_dic)))

    # create nodes
    revDef_dic = Dict(y => x[1] for x in collect(def_dic) for y in x[2])
    node_dic = Dict{Tuple{Vararg{Int}}, OptiNode}()
    foreach(x -> node_dic[x] = add_node(subGraph_dic[revDef_dic[x]]), plasmo_obj.nodes)

    # ! create variables
    # create dictionary for variable nameing
    strSets_dic = Dict{Symbol, Dict{Int, String}}()
    foreach(x -> strSets_dic[x] = Dict(v => string(k) for (k,v) in plasmo_obj.sets[x]), keys(plasmo_obj.sets))

    # initialize variable dictionary
    varSub_dic = plasmo_obj.varNode
    varDf_dic = plasmo_obj.var
    var_dic = Dict{Tuple{Int,Tuple{Vararg{Int}}}, Union{NodeVariableRef,VariableRef}}()

    for v in keys(varDf_dic)
        var_df = varDf_dic[v]
        # add subgraph info to variables
        var_df[!,:subGraph] .= map(x -> x in keys(varSub_dic) ? varSub_dic[x] : tuple(), var_df[!,:key])
        var_df = flatten(var_df, :subGraph)
        # add variables as a column
        setCol_arr = filter(x -> !(x in (:key,:upper,:lower,:subGraph)), namesSym(var_df))
        if wrtStr
            var_df[!,:var] = map(x -> createVar(string(v) * "[" *  join(map(y -> string(y) * ":" * strSets_dic[y][getindex(x,y)], setCol_arr),", ") * "]", x.lower, x.upper, isPlasmo ? node_dic[x.subGraph] : main_graph), eachrow(var_df))
        else
            var_df[!,:var] = map(x -> createVar(string(v) * "[" *  join(map(y -> getindex(x,y), setCol_arr),",") * "]", x.lower, x.upper, isPlasmo ? node_dic[x.subGraph] : main_graph), eachrow(var_df))
        end
        foreach(x -> var_dic[(x.key,x.subGraph)] = x.var, eachrow(var_df))
    end

    # ! create regular constraints
    cns_dic = plasmo_obj.cns
    for c in keys(cns_dic)
        cns_df = cns_dic[c]
        cns_df[!,:cnsEq] = map(x -> createCns(x.cnsObj, var_dic, x.subGraph, isPlasmo ? node_dic[x.subGraph] : main_graph), eachrow(cns_df))
    end

    # ! add linking constraints
    linkVar_df = DataFrame(node1 =Tuple[], node2 = Tuple[], var = Int[], cns = ConstraintRef[])
    for x in filter(x -> length(x[2]) > 1, collect(plasmo_obj.varNode))

        # get all involved subgraphs
        allSub_df = DataFrame(node = plasmo_obj.varNode[x[1]])
        allSub_df[!,:subGraph] .= map(y -> revDef_dic[y], allSub_df[!,:node])

        if length(unique(allSub_df[!,:subGraph])) > 1 # linking constraint across subgraphs

            # select top node within each subgraph for linking
            sort_df = combine(x -> (node = sortNodes(x.node)[1][1],), groupby(allSub_df, :subGraph))

            # sort overall nodes for linking
            sorted_arr = sortNodes(sort_df[!,:node])
            linkVar_df = createLinks(x[1], sorted_arr, linkVar_df, var_dic, main_graph)

        end

        # linking constraint within subgraph
        for z in groupby(allSub_df, :subGraph)
            relNodes_arr = collect(z.node)
            if length(relNodes_arr) > 1
                # determine structure of linking constraints based on tuples for nodes (= linking constraints should always be formulated with reference to a node higher in the hierarchy)
                groupNodes_arr = sortNodes(relNodes_arr)
                # create linking constraints looping over nodes
                linkVar_df = createLinks(x[1], groupNodes_arr, linkVar_df, var_dic, subGraph_dic[unique(allSub_df[!,:subGraph])[1]]; makeGra_boo = isPlasmo)
            end
        end

    end
    # ! write objective function
    # enfore objective on each node
    obj_df = DataFrame(node = Tuple{Vararg{Int}}[], obj = GenericAffExpr{Float64, NodeVariableRef}[])

    if isPlasmo
        for n in groupby(plasmo_obj.obj, :node)

            if plasmo_obj.objSense == :min
                obj_cns = @objective(node_dic[n.node[1]], Min, sum(map(x -> var_dic[(x.var, n.node[1])] * x.coeff, eachrow(n))))
            else
                obj_cns = @objective(node_dic[n.node[1]], Min, - sum(map(x -> var_dic[(x.var, n.node[1])] * x.coeff, eachrow(n))))
            end
            push!(obj_df, (node = n.node[1], obj = obj_cns))
        end
    else
        if plasmo_obj.objSense == :min
            @objective(main_graph, Min, sum(map(x -> var_dic[(x.var, x.node[1])] * x.coeff, eachrow(plasmo_obj.obj))))
        else
            @objective(main_graph, Min, -sum(map(x -> var_dic[(x.var, x.node[1])] * x.coeff, eachrow(plasmo_obj.obj))))
        end
    end
    # set graph objective accordingly
    foreach(x -> set_to_node_objectives(subGraph_dic[x]), keys(subGraph_dic))
    set_to_node_objectives(main_graph)

    # ! field declarations
    plasmo_obj.varMap = var_dic
    plasmo_obj.nodesMap = node_dic
    plasmo_obj.objFunc = obj_df
    plasmo_obj.linkCns = linkVar_df
    plasmo_obj.mainGraph = main_graph
    plasmo_obj.subGraphs = subGraph_dic

end

#endregion

#region # * handling of optimization problem

# ! create variable object
function createVar(name_str::String, lowBd_fl::Float64, upBd_fl::Float64, optModel::Union{Model,OptiNode}, bi::Bool = false)

	info = VariableInfo(!isnan(lowBd_fl), lowBd_fl, !isnan(upBd_fl), upBd_fl, false, NaN, false, NaN, bi, false)
	var_obj = JuMP.build_variable(error, info)
    JuMP.add_variable(optModel, var_obj, name_str)

end

# ! create actual constraint
function createCns(cns_obj::cnsObj, var_dic::Dict{Tuple{Int,Tuple{Vararg{Int}}}, Union{NodeVariableRef,VariableRef}}, subGraph_tup::Tuple{Vararg{Int}}, optModel::Union{Model,OptiNode})

    cns_expr = sum([var_dic[(x[1],subGraph_tup)] * x[2] for x in cns_obj.var])

    if cns_obj.sense == :(=)
        c_obj = @constraint(optModel, cns_expr == cns_obj.rhs)
    elseif cns_obj.sense == :(>=)
        c_obj = @constraint(optModel, cns_expr >= cns_obj.rhs)
    else cns_obj.sense == :(<=)
        c_obj = @constraint(optModel, cns_expr <= cns_obj.rhs)
    end

    return c_obj

end

# ! create linking constraints between variables of different nodes
function createLinks(var_int::Int, groupNodes_arr::Vector{<:Vector{<:Tuple{Vararg{Int64}}}}, linkVar_df::DataFrame, var_dic::Dict{Tuple{Int,Tuple{Vararg{Int}}}, Union{NodeVariableRef,VariableRef}}, graph_obj::Union{OptiGraph, Model}; makeGra_boo::Bool = true)

    # create linking constraints looping over nodes
    for g in eachindex(groupNodes_arr)
        # connect first variable on the level with previous level (if exists)
        if g != 1
            if makeGra_boo
                cns_obj = @linkconstraint(graph_obj, var_dic[(var_int, groupNodes_arr[g-1][1])] == var_dic[(var_int, groupNodes_arr[g][1])])
            else
                cns_obj = @constraint(graph_obj, var_dic[(var_int, groupNodes_arr[g-1][1])] == var_dic[(var_int, groupNodes_arr[g][1])])
            end
            push!(linkVar_df, (node1 = groupNodes_arr[g-1][1], node2 = groupNodes_arr[g][1], var = var_int, cns = cns_obj))
        end

        # connect first variable on the level with others on the same level
        if length(groupNodes_arr[g]) > 1
            for i in 2:length(groupNodes_arr[g])
                if makeGra_boo
                    cns_obj = @linkconstraint(graph_obj, var_dic[(var_int, groupNodes_arr[g][1])] == var_dic[(var_int, groupNodes_arr[g][i])])
                else
                    cns_obj = @constraint(graph_obj, var_dic[(var_int, groupNodes_arr[g][1])] == var_dic[(var_int, groupNodes_arr[g][i])])
                end
                push!(linkVar_df, (node1 = groupNodes_arr[g][1], node2 = groupNodes_arr[g][i], var = var_int, cns = cns_obj))
            end
        end

    end

    return linkVar_df
end

# ! detect conflicts in the model causing infeasibility (only available when not converted to Plasmo and using Gurobi)
function checkISS(plasmo_obj::lin2plasObj)

    if !isdefined(plasmo_obj, :mainGraph) || plasmo_obj.mainGraph isa OptiGraph
        error("Model was converted into Plasmo which does not support the IIS feature. Set isPlasmo=false in createOptProblem! to enable IIS.")
    end

    compute_conflict!(plasmo_obj.mainGraph)

    for cns in collect(plasmo_obj.cns)

        cns_arr = cns[2][!,:cnsEq]

        allConstr_arr = findall(map(x -> MOI.ConflictParticipationStatusCode(0) != MOI.get(plasmo_obj.mainGraph.moi_backend, MOI.ConstraintConflictStatus(), x.index), cns_arr))

        if !isempty(allConstr_arr)
            println("$(length(allConstr_arr)) of IIS in $(cns[1]) constraints.")
            for iisConstr in allConstr_arr
                println(cns[2][iisConstr,:])
            end
        end

    end

end

#endregion

#region # * low-level data processing

# ! convert python array to julia array of provided data type
convArray(inArr_py::Py, dataType::Type) = map(x -> pyconvert(dataType,x), inArr_py)
convXArray(inXArr_py::Py, dataType::Type) = map(x -> pyconvert(dataType,x), inXArr_py.flatten())
convXArray(inXArr_py::Py, dataType::Type, reIdx_arr::Vector{Int}) = map(x -> pyconvert(dataType,x), inXArr_py.flatten()[reIdx_arr])

namesSym(in_df::DataFrame) = Symbol.(names(in_df))

# ! convert constraint xarray to dataframe
function convertXarr(filt_xarr::Py, value_arr::Vector{String}, label_arr::Vector{String}, sign_boo::Bool = false)

    stack_xarr = filt_xarr.stack(z=filt_xarr.dims).dropna("z")

    # convert set data to dictionary
    setData_dic = Dict{Symbol, Vector{String}}()
    for set in filter(x -> occursin("set", string(x)), string.(collect(stack_xarr.coords.keys())))
        setData_dic[Symbol(set)] = string.(pyconvert(Vector, stack_xarr[set].values))
    end

    # convert value data to dictionary
    setValue_dic = Dict{Symbol, Vector{Float64}}()
    for var in value_arr
        setValue_dic[Symbol(var)] = pyconvert(Vector, stack_xarr[var].values)
    end

    # convert label data to dictionary
    setLabel_dic = Dict{Symbol, Vector{Int}}()
    for label in label_arr
        setLabel_dic[Symbol(label)] = Integer.(pyconvert(Vector, stack_xarr[label].values))
    end

    # write to a dataframe
    dataCns_df = DataFrame(setData_dic)
    foreach(x -> dataCns_df[!, x] = setValue_dic[x], keys(setValue_dic))
    foreach(x -> dataCns_df[!, x] = setLabel_dic[x], keys(setLabel_dic))
    if sign_boo dataCns_df[!,:sign] = string.(pyconvert(Vector, stack_xarr["sign"].values)) end

    return dataCns_df
end

# ! update sets_dic with new values
function updateSetsDic!(sets_dic::Dict{Symbol, Dict{String, Int}}, setVal_dic::Dict{Symbol, Vector{String}})

    for s in keys(setVal_dic)
        if s in keys(sets_dic) # add new data to existing sub-dictonary
            sub_dic = sets_dic[s]
            len_int = length(sub_dic)
            foreach(y -> sub_dic[setVal_dic[s][y]] = len_int + y, eachindex(collect(setdiff(setVal_dic[s], keys(sub_dic)))))
        else # create new sub
            sets_dic[s] = Dict(setVal_dic[s][x] => x for x in eachindex(setVal_dic[s]))
        end
    end

end

# ! replace set columns of strings with integer indices
function replaceSetColumns(in_df::DataFrame, exSet_dic::Dict{Symbol, Vector{String}}, sets_dic::Dict{Symbol, Dict{String, Int}})
    newCol_arr = Symbol.(keys(exSet_dic),:_)
    in_df = rename(in_df, Symbol.(keys(exSet_dic)) .=> newCol_arr)
    foreach(s -> in_df[:,s] = map(y -> sets_dic[s][typeof(y) == String ? y : string(y)], in_df[:,Symbol(s,:_)]), collect(keys(exSet_dic)))
    select!(in_df, Not(newCol_arr))
end

# ! replace set columns of integer indices with strings
function replaceSetColumns(in_df::DataFrame, revSets_dic::Dict{Symbol, Dict{Int, String}})
    exCol_arr = intersect(namesSym(in_df), keys(revSets_dic))
    newCol_arr = Symbol.(exCol_arr,:_)
    in_df = rename(in_df, exCol_arr .=> newCol_arr)
    foreach(s -> in_df[:,s] = map(y -> revSets_dic[s][y], in_df[:,Symbol(s,:_)]), exCol_arr)
    select!(in_df, Not(newCol_arr))
end

# ! sort nodes according to hierarchy
function sortNodes(in_arr::Union{Array{<:Tuple{Vararg{Int64}}},SubArray{<:Tuple{Vararg{Int64}}}})
    sortNodes_arr = map(x -> (x,findall(x .== 0)), sort(in_arr, rev = true)) |> (y -> sort(y, by = x -> x[2], rev = true))
    out_arr = [getindex.(filter(y -> y[2] == x, sortNodes_arr),1) for x in unique(getindex.(sortNodes_arr,2))]
    return out_arr
end

# ! write infos on nodes into dictionary of variables
function addToVarSubDic!(varSub_dic::Dict{Int, Array{Tuple{Vararg{Int}}}}, var_set::Vector{T}, subGraph_tup::Tuple{Vararg{Int}}) where {T <: Integer}

    for v in var_set
        if v in keys(varSub_dic) # add to existing entry
            if !(subGraph_tup in varSub_dic[v])
                push!(varSub_dic[v], subGraph_tup)
            end
        else # create new entry
            varSub_dic[v] = [subGraph_tup]
        end
    end

end

#endregion
