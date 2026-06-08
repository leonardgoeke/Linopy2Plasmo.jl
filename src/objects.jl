mutable struct cnsObj
	var::Vector{Pair{Integer, Float64}}
    rhs::Float64
    sense::Symbol
    function cnsObj(var::Vector{Pair{T, Float64}}, rhs::Float64, sense::Symbol) where {T <: Integer}
        new(var, rhs, sense)
    end
end

mutable struct lin2plasObj
	# problem structure
    sets::Dict{Symbol, Dict{String, Int}} # dictionary of sets
    revSets::Dict{Symbol, Dict{Int, String}} # dictionary of sets with reversed keys and values
    setHier::Tuple{Vararg{Symbol}} # hierarchy of sets for defining nodes
    nodes::Array{Tuple{Vararg{Int}}} # array of all nodes
    nodesMap::Dict{Tuple{Vararg{Integer}}, OptiNode} # maps nodes to OptiNode objects
    mainGraph::OptiGraph # main graph of the problem
    subGraphs::Dict{Symbol, OptiGraph} # dictionary of subgraphs
    # variables
    var::Dict{Symbol, DataFrame} # all variables stored in dataframe
    varNode::Dict{Int,Array{Tuple{Vararg{Int}}}} # maps each variables to relevant nodes
    varMap::Dict{Tuple{Int,Tuple{Vararg{Int}}}, Union{NodeVariableRef,VariableRef}} # maps variable and node id to variable
    linkCns::DataFrame # linking constraints between nodes and subgraphs
    # constraints and objective
    cns::Dict{Symbol, DataFrame} # all constraints stored in dataframe
    obj::DataFrame # objective function stored in dataframe
    objSense::Symbol # objective function sense
    objFunc::DataFrame # objectives of each node

    lin2plasObj() = new()
    function lin2plasObj(linopyMod_str::String)

        plasmo_obj = new()

        # ! get model file
        model = linopy.io.read_netcdf(linopyMod_str)
        println("Read-in .netcdf file")

        # ! read-in sets and variables
        sets_dic = Dict{Symbol, Dict{String, Int}}()
        varDf_dic = Dict{Symbol, DataFrame}()

        println("Loop over variables")
        start = now()
        for v in [string(k) for (k, _) in model.variables.items()]

            # store variabeles as dataframe
            data_py = model.variables[v].data
            filt_xarr = data_py.where((data_py.labels != -1), drop = true)
            var_df = convertXarr(filt_xarr, ["lower", "upper"], ["labels"])

            # update set dictonary
            exSet_dic = Dict(x => unique(var_df[!,x]) |> (z -> typeof(z) >: Vector{String} ? z : string.(z)) for x in filter(y -> !(y in (:lower,:upper,:labels)), namesSym(var_df)))
            updateSetsDic!(sets_dic, exSet_dic)

            # replace set columns of strings with integer indices
            var_df = replaceSetColumns(var_df, exSet_dic, sets_dic)

            varDf_dic[Symbol(v)] = rename(var_df,:labels => :key)

        end
        println("Variables processed in $(now() - start)")

        # ! read-in constraints
        cns_dic = Dict{Symbol, DataFrame}()

        println("Loop over constraints")
        start = now()
        for c in [string(k) for (k, _) in model.constraints.items()]

            # store constraints as dataframe
            data_py = model.constraints[c].data
            filt_xarr = data_py.where((data_py.vars != -1) .& (data_py.coeffs != 0.0) .& (data_py.rhs != Inf), drop = true)
            cns_df = convertXarr(filt_xarr, ["coeffs", "rhs"], ["vars", "labels"] , true)
            if isempty(cns_df) continue end

            # update set dictonary
            exSet_dic = Dict(x => unique(cns_df[!,x]) |> (z -> typeof(z) >: Vector{String} ? z : string.(z)) for x in filter(y -> !(y in (:_term,:coeffs,:vars,:sign,:rhs,:labels)), namesSym(cns_df)))
            updateSetsDic!(sets_dic, exSet_dic)

            # replace set columns of strings with integer indices
            cns_df = replaceSetColumns(cns_df, exSet_dic, sets_dic)

            # group to get constraint again
            set_arr = collect(keys(exSet_dic))
            setCns_arr = vcat(set_arr,[:cnsObj])
            cns_df = combine(x -> (; zip(setCns_arr, vcat(map(y -> x[1,y],set_arr),[cnsObj((x.vars .=> x.coeffs), x.rhs[1], Symbol(x.sign[1]))]))...), groupby(cns_df, :labels))

            cns_dic[Symbol(c)] = cns_df

        end
        println("Constraints processed in $(now() - start)")

        # ! extract info on objective function
        # get variables of objective function and to which node version to enforce them
        var_arr = map(x -> pyconvert(Int,x), model.objective.expression.vars.data)

        # get coefficients of objective function
        coeffs_arr = map(x -> pyconvert(Float64,x), model.objective.expression.coeffs.data)
        obj_df = DataFrame(var = getindex.(var_arr,1), coeff = coeffs_arr)

        plasmo_obj.sets = sets_dic
        plasmo_obj.revSets = Dict(s => Dict(x[2] => x[1] for x in collect(sets_dic[s])) for s in keys(sets_dic))
        plasmo_obj.var = varDf_dic
        plasmo_obj.cns = cns_dic
        plasmo_obj.obj = obj_df
        plasmo_obj.objSense = Symbol(model.objective.sense)


        return plasmo_obj
    end
end
