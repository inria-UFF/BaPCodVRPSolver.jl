module RCSP

# package code goes here
using JuMP
using LightGraphs

using MathProgBase.SolverInterface
using ..BlockDecompositionExtras.BlockSolverExtraInterface
import ..BlockDecompositionExtras.BlockSolverExtraInterface:
  send_rcsp_to_solver!

import LightGraphs.add_edge!
import LightGraphs.vertices

# types
export Network

# user methods
export  add_edge!,
        set_edge_vars!,
        addresource!,
        addmainresource!,
        addspecialresource!,
        vertexconsumptionbounds!,
        vertexspecialconsumptionbounds!,
        edgeconsumption!,
        arcconsumptionbounds!,
        edgespecialconsumption!,
        attach_elem_set_to_vertex!,
        attach_elem_set_to_edge!,
        add_vertex_to_mem_of_elem_set!,
        add_edge_to_mem_of_elem_set!,
        add_vertex_to_packing_set!,
        add_edge_to_packing_set!,
        add_vertex_to_covering_set!,
        add_edge_to_covering_set!,
        add_packing_set_to_packing_set_cut_neighbourhood!,
	    set_elementarity_sets_distance_matrix!,
        add_rcsp_capacity_cuts!,
        add_rcsp_strongkpath_cuts!,
        add_rcsp_clique_cuts!,
        add_rcsp_lim_mem_rank1_cuts!,
        generate_rcsp_oracle_for_dw_sp!,
        vertices,
        pos2edge, edge2pos,
        add_rcsp_paths_per_network_branching!,
        add_rcsp_elemset_assign_branching!,
        add_rcsp_elemset_resource_consumption_branching!,
        add_rcsp_packset_ryanfoster_branching!,
        add_rcsp_knapsack_consumption_cut!,
        associate_var_to_resource!,
        save_standalone!,
        add_permanent_ryan_foster_constraint!

# Send data to BaPCod
function send_rcsp_to_solver!(model)
  set_rcsp_oracles!(model.solver, model.ext[:RCSP_oracles])
  if haskey(model.ext, :RCSP_generic_cuts)
    set_rcsp_cuts!(model.solver, model.ext[:RCSP_generic_cuts])
  end
  if haskey(model.ext, :RCSP_branching)
    set_rcsp_branching!(model.solver, model.ext[:RCSP_branching])
  end
end

mutable struct Network
  nb_nodes::Integer
  digraph
  source::Integer
  sink::Integer
  nb_elementaritysets::Integer
  nb_packingsets::Integer
  nb_coveringsets::Integer
  #elemsets_id::Vector{Integer} # es_pos => es_id
  elemsets_ni_nodes::Vector{Vector{Integer}} # es_pos to nodes
  elemsets_ni_edges::Vector{Vector{Integer}} # es_pos to edges
  packsets_ni_nodes::Vector{Vector{Integer}} # ps_pos => nodes
  packsets_ni_edges::Vector{Vector{Integer}}
  covsets_ni_nodes::Vector{Vector{Integer}}
  covsets_ni_edges::Vector{Vector{Integer}}
  moes_vertex::Vector{Vector{Integer}} # es_pos => nodes
  moes_edges::Vector{Vector{Integer}} # es_pos => edges
  ps2pscutneighbourhood::Vector{Vector{Integer}} # ps => ps2add
  edge_id2pos::Dict{Tuple, Integer} # (s, d, o) => edge_pos
  edge_pos2id::Vector{Tuple} # edge_pos => (s, d, o)
  bcid_2pos::Dict{Int,Int}
  multiedges::Dict{Tuple, Integer} # (s, d) => o_max
  edgescosts::Vector{Cdouble} # edge_pos => cost
  edgesvars::Vector{Vector{Tuple}} # edge_pos => variable
  resources::Dict{Integer, Tuple}
  specialresources::Vector{Bool}
  verticesconsumptionsbounds::Array{Dict{Integer, Tuple}} # [n][r] => (lb, ub)
  arcsconsumptionsbounds::Array{Dict{Integer, Tuple}} # [a][r] => (lb, ub)
  verticesspecialconsumptionsbounds::Array{Dict{Integer, Tuple}} # [n][r] => (lb, ub)
  edgesconsumptions::Vector{Dict{Integer, Cdouble}} #[edge_pos][r] => val
  edgesspecialconsumptions::Vector{Dict{Integer, Cdouble}} #[edge_pos][r] => val
  resource_associated_var::Dict{Integer, Integer} # r => var
  save_standalone::Bool
  standalone_filename::String
  c_net
  es_dist_matrix
  ryanfoster_constraints::Vector{Tuple{Integer,Integer,Bool}} # (firstPackSetId,secondPackSetId,together)
end

function doublevector_generator(size)
  return [Vector{Int}() for i in 1:size]
end
"""
    Network(n::Integer)

Create a Network with `n` nodes. The set of vertices is ``⟦1,n⟧``.
The source and the sink node are respectively nodes `1` and `n`.


**Optional arguments**

    Network(n::Integer, source = s, sink = d, elemsets = e, packsets = p, covsets = c)

Nodes `s` and `d` are respectively the source and the sink.
The number of elementarity sets is `e`, the number of packing sets is `p` and
the number of covering sets is `c`.

These values are zero by default.
"""
function Network(nb_nodes::Integer; source = 1, sink = nb_nodes, elemsets = 0, packsets = 0, covsets = 0)
  # Check the source & the sink
  (source < 1 || source > nb_nodes) && error("Unknown source node $source.")
  (sink < 1 || sink > nb_nodes) && error("Unknown sink node $sink.")
  # Init vertices consumption array{array}
  init_vc = Array{Dict}(undef, nb_nodes)
  for i in 1:length(init_vc)
    init_vc[i] = Dict{Integer, Tuple}()
  end
  # Create the object
  return Network(nb_nodes,
                 DiGraph(nb_nodes),
                 source, sink, elemsets, packsets, covsets,
                 doublevector_generator(elemsets),
                 doublevector_generator(elemsets),
                 doublevector_generator(packsets),
                 doublevector_generator(packsets),
                 doublevector_generator(covsets),
                 doublevector_generator(covsets),
                 doublevector_generator(elemsets), # should be replaced by a dict
                 doublevector_generator(elemsets), # should be replaced by a dict
                 doublevector_generator(packsets), # should be replaced by a dict
                 Dict{Tuple, Integer}(),
                 Vector{Tuple}(),
		 Dict{Int,Int}(),
                 Dict{Tuple, Integer}(),
                 Vector{Cdouble}(),
                 Vector{Vector}(),
                 Dict{Integer, Tuple}(),
                 Bool[],
                 init_vc,
		 deepcopy(init_vc),
		 deepcopy(init_vc),
                 Vector{Dict}(),
                 Vector{Dict}(),
                 Dict{Integer, Integer}(),
                 false,
                 "empty",
                 nothing,
		 nothing, 
     Vector{Tuple{Integer,Integer,Bool}}()
                 )
end

function check_node(net::Network, node::Integer)
  (node > net.nb_nodes) && error("Unknown node $node.")
end

function check_edge_pos(net::Network, edge_pos::Integer)
  (edge_pos > length(net.edge_id2pos)) && error("Unknown edge position $edge_pos.")
end

# Translate functions
function pos2edge(net::Network, edge_pos::Integer)
  check_edge_pos(net, edge_pos)
  return net.edge_pos2id[edge_pos]
end

edge2pos(net::Network, edge::Tuple) = net.edge_id2pos[edge]

#auxiliary function used by add_edge! and set_edge_vars!
function check_var_info(varinfo)
    if typeof(varinfo) <: Tuple
      if typeof(varinfo[1]) != JuMP.Variable
        error("The first element of a Tuple must be a JuMP.Variable.")
      end
      if !(typeof(varinfo[2]) <: Real)
        error("The second element of a Tuple must be a Real.")
      end
    elseif typeof(varinfo) != JuMP.Variable
      error("Elements of the variable array must be Tuple or JuMP variables.")
    end
end


"""
    add_edge!(net::Network, queue::Integer, tail::Integer)

Add the arc ``(queue, tail)`` to the network `net` and returns an integer which
is the id of the arc.

    add_edge!(net::Network, arc::Tuple)

Add the arc `arc` to the network `net` and returns an integer which is the id of
the arc.

**Optional arguments**

    add_edge!(net::Network, queue::Integer, tail::Integer, cost = c, var = x)

The initial cost of the arc is `c` and the JuMP variable attached to the arc is `x`.
In the pricing subproblem, the cost
of the arc is the sum of the initial cost and the cost of the variable.

By default, the initial cost of an arc is ``0`` and no variable is attached.
"""
function add_edge!(net::Network, s::Integer, d::Integer; cost = 0, var = nothing)
  # We add the arc if it does not exist
  if !haskey(net.multiedges, (s, d))
    status = add_edge!(net.digraph, s, d)
    !status && error("Cannot create edge between $s and $d.")
    net.multiedges[(s, d)] = 0
  end
  # In case of multiple arcs
  net.multiedges[(s, d)] += 1
  o = net.multiedges[(s, d)]
  # Compute the id of the arc
  edge_pos = length(net.edge_id2pos) + 1
  net.edge_id2pos[(s, d, o)] = edge_pos
  push!(net.edge_pos2id, (s, d, o))
  # Cost of the arc
  push!(net.edgescosts, Cdouble(cost))

  varsinfo = Vector{Tuple}() # array of (column, coeff)
  if var != nothing
    if typeof(var) == JuMP.Variable
      push!(varsinfo, (var.col, 1))
    elseif typeof(var) <: Tuple
      check_var_info(var)
      push!(varsinfo, (var[1].col, var[2]))
    elseif typeof(var) <: Array
      for v in var
        check_var_info(v)
        if typeof(v) == JuMP.Variable
          push!(varsinfo, (v.col, 1))
        else
          push!(varsinfo, (v[1].col, v[2]))
        end
      end
    else
      @warn("Type of the variable attached to the arc ($s, $d) is $(typeof(var)). Only JuMP.Variable, Tuple or Array are supported.")
    end
  end
  push!(net.edgesvars, varsinfo)
  # Init arc consumptions
  push!(net.edgesconsumptions, Dict{Integer, Cdouble}())
  push!(net.edgesspecialconsumptions, Dict{Integer, Cdouble}())
  push!(net.arcsconsumptionsbounds, Dict{Integer, Tuple}())
  return edge_pos
end

function add_edge!(net::Network, edge::Tuple; cost = 0, var = nothing)
  s, d = edge
  return add_edge!(net, s, d, cost = cost, var = var)
end

"""
    set_edge_vars!(net::Network, edge_pos::Integer, vars)

Sets the variables associated to the edge
"""
function set_edge_vars!(net::Network, edge_pos, vars)

  check_edge_pos(net, edge_pos)

  varsinfo = Vector{Tuple}() # array of (column, coeff)
  if vars != nothing
    if typeof(vars) == JuMP.Variable
      push!(varsinfo, (vars.col, 1))
    elseif typeof(vars) <: Tuple
      check_var_info(vars)
      push!(varsinfo, (vars[1].col, var[2]))
    elseif typeof(vars) <: Array
      for v in vars
        check_var_info(v)
        if typeof(v) == JuMP.Variable
          push!(varsinfo, (v.col, 1))
        else
          push!(varsinfo, (v[1].col, v[2]))
        end
      end
    else
      @warn("Type of the variable attached to the arc is $(typeof(vars)). Only JuMP.Variable, Tuple or Array are supported.")
    end
  end
  net.edgesvars[edge_pos] = varsinfo
end


"""
    addresource!(net::Network, r::Integer)

Create a resource with the id `r` in the network `net`.
"""
function addresource!(net::Network, res_id::Integer; disposable = true)
  if !haskey(net.resources, res_id)
    net.resources[res_id] = (false, disposable, 0)
  else
    error("Resource with id $res_id already exists.")
  end
  return res_id
end

"""
    addmainresource!(net::Network, r::Integer, stepsize = ss)

Create a main resource with the id `r` in the network `net`.
The step size argument is optional and its default value is `1.0`.
"""
function addmainresource!(net::Network, res_id::Integer; stepsize=1.0, disposable = true)
  if !haskey(net.resources, res_id)
    net.resources[res_id] = (true, disposable, Cdouble(stepsize))
  else
    error("Resource with id $res_id already exists.")
  end
  return res_id
end

function addspecialresource!(net::Network; disposable = true)
  if length(net.specialresources) < 512
    push!(net.specialresources, disposable)
  else
    error("Maximum number of special resource is 512.")
  end
  return length(net.specialresources)
end


### Vertex resources consumptions
"""
    vertexconsumption!(net::Network, v::Integer, r::Integer; lb = l, ub = u, value = c)

Set the consumption of resource `r` at vertex `v` in network `net`.

*Optional arguments* :
 - Parameter `value` defines the resource consumption of the vertex.
 - Parameters `lb` and `ub` define lower bound and the upper bound of resource consumption at the vertex.
"""
function vertexconsumptionbounds!(net::Network, node_id::Integer, res_id::Integer; kwargs...)
  (node_id > net.nb_nodes) && error("Unknown vertex $node_id.")
  throwerror = true
  lb, ub = NaN, NaN
  dict = Dict(k => v for (k, v) in kwargs)
  if haskey(dict, :lb)
    lb = Cdouble(dict[:lb])
    throwerror = false
  end
  if haskey(dict, :ub)
    ub = Cdouble(dict[:ub])
    throwerror = false
  end
  net.verticesconsumptionsbounds[node_id][res_id] = (lb, ub)
  throwerror && error("vertexconsumption! : invalid parameters. Use `lb`or `ub`.")
end

function vertexspecialconsumptionbounds!(net::Network, node_id::Integer, res_id::Integer; kwargs...)
  (node_id > net.nb_nodes) && error("Unknown vertex $node_id.")
  throwerror = true
  lb, ub = NaN, NaN
  dict = Dict(k => v for (k, v) in kwargs)
  if haskey(dict, :lb)
    lb = Cdouble(dict[:lb])
    throwerror = false
  end
  if haskey(dict, :ub)
    ub = Cdouble(dict[:ub])
    throwerror = false
  end
  net.verticesspecialconsumptionsbounds[node_id][res_id] = (lb, ub)
  throwerror && error("vertexspecialconsumption! : invalid parameters. Use `lb`or `ub`.")
end


### Edge resources consumptions and bounds
"""
    edgeconsumption!(net::Network, a::Integer, r::Integer; ; lb = l, ub = u, value = val)

Set the amount of resource `r` consumed by the arc `a` in network `net`.

*Optional arguments* :
 - Parameter `value` defines the resource consumption of the arc.
 - Parameters `lb` and `ub` define lower bound and the upper bound of resource consumption at the arc.
"""
function edgeconsumption!(net::Network, edge_pos::Integer, res_id::Integer; kwargs...)
  dict = Dict(k => v for (k, v) in kwargs)
  if haskey(dict, :value)
    net.edgesconsumptions[edge_pos][res_id] = Cdouble(dict[:value])
  else
    error("edgeconsumption! only takes the parameter `value`.")
  end
end

function edgeconsumption!(net::Network, edge::Tuple, res_id::Integer; kwargs...)
  edgeconsumption!(net, edge2pos(edge), res_id, kwargs)
end

function edgespecialconsumption!(net::Network, edge_pos::Integer, res_id::Integer; kwargs...)
  dict = Dict(k => v for (k, v) in kwargs)
  if haskey(dict, :value)
    net.edgesspecialconsumptions[edge_pos][res_id] = Cdouble(dict[:value])
  else
    error("edgeconsumption! only takes the parameter `value`.")
  end
end

function edgespecialconsumption!(net::Network, edge::Tuple, res_id::Integer; kwargs...)
  edgespecialconsumption!(net, edge2pos(edge), res_id, kwargs)
end

function arcconsumptionbounds!(net::Network, arc_id::Integer, res_id::Integer; kwargs...)
  throwerror = true
  lb, ub = NaN, NaN
  dict = Dict(k => v for (k, v) in kwargs)
  if haskey(dict, :lb)
    lb = Cdouble(dict[:lb])
    throwerror = false
  end
  if haskey(dict, :ub)
    ub = Cdouble(dict[:ub])
    throwerror = false
  end
  net.arcsconsumptionsbounds[arc_id][res_id] = (lb, ub)
  throwerror && error("arcconsumptionbounds! : invalid parameters. Use `lb`or `ub`.")
end

function check_elemset(net::Network, es_id)
  (0 < es_id <= net.nb_elementaritysets) || error("Unknown elementarity set $es_id.")
end

function check_packset(net::Network, ps_id)
  (0 < ps_id <= net.nb_packingsets) || error("Unknown packing set $ps_id.")
end

function check_covset(net::Network, cs_id)
  (0 < cs_id <= net.nb_coveringsets) || error("Unknown covering set $cs_id.")
end

"""
    attach_elem_set_to_vertex!(net::Network, v::Integer, es::Integer)

Attach the elementarity set `es` to the vertex `v` in the network `net`.
"""
function attach_elem_set_to_vertex!(net::Network, n_id::Integer, es_id::Integer)
  (n_id > net.nb_nodes) && error("Unknown node $n_id.")
  check_elemset(net, es_id)
  push!(net.elemsets_ni_nodes[es_id], n_id)
end

"""
    attach_elem_set_to_edge!(net::Network, a::Integer, es::Integer)

Attach the elementarity set `es` to the arc `a` in the network `net`.
"""
function attach_elem_set_to_edge!(net::Network, edge_pos::Integer, es_id::Integer)
  check_edge_pos(net, edge_pos)
  check_elemset(net, es_id)
  push!(net.elemsets_ni_edges[es_id], edge_pos)
end


"""
    add_vertex_to_mem_of_elem_set!(net::Network, v::Int, es::Int)

Add the vertex `v` to the memory of the elementarity set `es` in network `net`.
"""
function add_vertex_to_mem_of_elem_set!(net::Network, n_id::Int, es_id::Int)
  (n_id > net.nb_nodes) && error("Unknown node $n_id")
  check_elemset(net, es_id)
  push!(net.moes_vertex[es_id], n_id)
end

"""
    add_edge_to_mem_of_elem_set!(net::Network, edge_pos::Int, es_id::Int)

Add the edge with the id `edge_pos` the memory of the elementarity set `es_id`.
"""
function add_edge_to_mem_of_elem_set!(net::Network, edge_pos::Int, es_id::Int)
  check_elemset(net, es_id)
  push!(net.moes_edges[es_id], edge_pos)
end

"""
    add_vertex_to_packing_set!(net::Network, v::Int, ps::Int)

Add the vertex `v` to the packing set `ps` in network `net`.
"""
function add_vertex_to_packing_set!(net::Network, n_id::Int, ps_id::Int)
  check_packset(net, ps_id)
  push!(net.packsets_ni_nodes[ps_id], n_id)
end

"""
    add_edge_to_packing_set!(net::Network, a::Int, ps::Int)

Add arc `a` to the packing set `ps` in network `net`.
"""
function add_edge_to_packing_set!(net::Network, edge_pos::Int, ps_id::Int)
  check_packset(net, ps_id)
  push!(net.packsets_ni_edges[ps_id], edge_pos)
end

"""
    add_vertex_to_covering_set!(net::Network, v::Int, cs::Int)

Add vertex `v` to the covering set `cs` in network `net`.
"""
function add_vertex_to_covering_set!(net::Network, n_id::Int, cs_id::Int)
  check_covset(net, cs_id)
  push!(net.covsets_ni_nodes[cs_id], n_id)
end

"""
    add_edge_to_covering_set!(net::Network, a::Integer, cs::Int)

Add arc `a` to the covering set `cs` in network `net`.
"""
function add_edge_to_covering_set!(net::Network, edge_pos, cs_id::Int)
  check_covset(net, cs_id)
  push!(net.covsets_ni_edges[cs_id], edge_pos)
end

"""
    add_packing_set_to_packing_set_cut_neighbourhood!(net::Network, ps::Int, mps::Int)

Add packing set `ps` to the memory of packing set `mps` in network `net`.
"""
function add_packing_set_to_packing_set_cut_neighbourhood!(net::Network, ps_id_to_add::Int, ps_id::Int)
  check_packset(net, ps_id_to_add)
  check_packset(net, ps_id)
  push!(net.ps2pscutneighbourhood[ps_id], ps_id_to_add)
end

function set_elementarity_sets_distance_matrix!(net::Network, matrix::Array{Array{Float64,1},1})
  nb_es = net.nb_elementaritysets
  if length(matrix) != nb_es
     error("distance matrix size does not equal to the number of elementarity sets")
     for row in matrix
        if length(row) != nb_es
           error("distance matrix size does not equal to the number of elementarity sets")
	end
     end
  end
  cdouble_matrix = [[Cdouble(i) for i in row] for row in matrix]
  net.es_dist_matrix = cdouble_matrix
end

function add_permanent_ryan_foster_constraint!(net::Network, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool)
  push!(net.ryanfoster_constraints, (firstPackSetId, secondPackSetId, together))
end

# """
#     add_elem_set_to_mem_of_elem_set_cut!(net::Network, es_id_2add::Int, es_id::Int)
#
# Add the elementarity set `es_id_2add` to the memory of the elementarity set `es_id`.
# """
# function add_elem_set_to_mem_of_elem_set_cut!(net::Network, es_id_2add::Int, es_id::Int)
#   id = findfirst(net.elemsets_id, es_id)
#   (id == 0) && error("Elementarity set with the id $es_id does not exist.")
#   id_2add = findfirst(net.elemsets_id, es_id_2add)
#   (id_2add == 0) && error("Elementarity set with the id $es_id_2add does not exist.")
#   push!(net.moes_elemset[id], id_2add)
# end

### Cuts
"""
    add_rcsp_capacity_cuts!(model::JuMP.Model, max_capacity, demands,
                          is_facultative = true, root_priority_level = 1.0,
                          non_root_priority_level = 1.0)

Add capacity cut
"""
function add_rcsp_capacity_cuts!(model::JuMP.Model, max_capacity, demands;
                            is_facultative = true, root_priority_level = 1.0,
                            non_root_priority_level = 1.0, two_path_cuts_res_id = -1)
  if !haskey(model.ext, :RCSP_generic_cuts)
    model.ext[:RCSP_generic_cuts] = Vector{Tuple}() # Pair (:type, args)
  end
  cint_demands = [Cint(d) for d in demands]
  args = (max_capacity, cint_demands, is_facultative, root_priority_level, non_root_priority_level, two_path_cuts_res_id)
  push!(model.ext[:RCSP_generic_cuts], (:capacity, args))
end

"""
    add_rcsp_strongkpath_cuts!(model::JuMP.Model, max_capacity, demands,
                               is_facultative = true, root_priority_level = 1.0,
                               non_root_priority_level = 1.0)

Add strong k-path cut
"""
function add_rcsp_strongkpath_cuts!(model::JuMP.Model, max_capacity, demands;
                                    is_facultative = true, root_priority_level = 1.0,
                                    non_root_priority_level = 1.0)
  if !haskey(model.ext, :RCSP_generic_cuts)
    model.ext[:RCSP_generic_cuts] = Vector{Tuple}() # Pair (:type, args)
  end
  cint_demands = [Cint(d) for d in demands]
  args = (max_capacity, cint_demands, is_facultative, root_priority_level, non_root_priority_level)
  push!(model.ext[:RCSP_generic_cuts], (:strongkpath, args))
end

"""
    add_rcsp_clique_cuts!(model::JuMP.Model)

Add clique cuts.
"""
function add_rcsp_clique_cuts!(model::JuMP.Model)
  if !haskey(model.ext, :RCSP_generic_cuts)
    model.ext[:RCSP_generic_cuts] = Vector{Tuple}() # Pair (:type, args)
  end
  push!(model.ext[:RCSP_generic_cuts], (:clique, nothing))
end

"""
    add_rcsp_lim_mem_rank1_cuts!(model::JuMP.Model)

Add limited memory rank 1 cuts.
"""
function add_rcsp_lim_mem_rank1_cuts!(model::JuMP.Model)
  if !haskey(model.ext, :RCSP_generic_cuts)
    model.ext[:RCSP_generic_cuts] = Vector{Tuple}() # Pair (:type, args)
  end
  push!(model.ext[:RCSP_generic_cuts], (:limmemrank1, nothing))
end

### register the RCSP oracle
"""
    generate_rcsp_oracle_for_dw_sp!(model::JuMP.Model, sp::Tuple, gen_rcsp_fct::Function)

Assign the oracle `gen_rcsp_fct` to the subproblem with the id `sp`.
The function `gen_rcsp_fct` is written by the user and returns a `Network`.
"""
function generate_rcsp_oracle_for_dw_sp!(model::JuMP.Model, sp::Tuple, gen_rcsp_fct::Function)
  # attach the function and the id of the dw sp to the model
  if !haskey(model.ext, :RCSP_oracles)
    model.ext[:RCSP_oracles] = Dict{Tuple, Function}()
  end
  (sp_type, sp_id) = sp
  if isa(sp_id, Integer) sp_id = (sp_id,) end
  model.ext[:RCSP_oracles][(sp_type, sp_id)] = gen_rcsp_fct
end

function add_rcsp_paths_per_network_branching!(model::JuMP.Model, priority::Real)
  if !haskey(model.ext, :RCSP_branching)
    model.ext[:RCSP_branching] = Dict{Symbol, Any}()
  end
  model.ext[:RCSP_branching][:paths_per_network_branching] = priority
end

function add_rcsp_elemset_assign_branching!(model::JuMP.Model, priority::Real)
  if !haskey(model.ext, :RCSP_branching)
    model.ext[:RCSP_branching] = Dict{Symbol, Any}()
  end
  model.ext[:RCSP_branching][:elemset_assign_branching] = priority
end

function add_rcsp_elemset_resource_consumption_branching!(model::JuMP.Model, priority::Real)
  if !haskey(model.ext, :RCSP_branching)
    model.ext[:RCSP_branching] = Dict{Symbol, Any}()
  end
  model.ext[:RCSP_branching][:elemset_rescons_branching] = priority
end

function add_rcsp_packset_ryanfoster_branching!(model::JuMP.Model, priority::Real)
  if !haskey(model.ext, :RCSP_branching)
    model.ext[:RCSP_branching] = Dict{Symbol, Any}()
  end
  model.ext[:RCSP_branching][:packset_ryanfoster_branching] = priority
end

function add_rcsp_knapsack_consumption_cut!(model::JuMP.Model, rootPrio::Real, nonRootPrio::Real)
  # TODO create a new map for cuts (or merge cuts & braching)
  if !haskey(model.ext, :RCSP_branching)
    model.ext[:RCSP_branching] = Dict{Symbol, Any}()
  end
  model.ext[:RCSP_branching][:knpconsumption_branching] = (rootPrio, nonRootPrio)
end

function associate_var_to_resource!(net::Network, res_id::Integer, var::JuMP.Variable)
  net.resource_associated_var[res_id] = var.col - 1
end

function save_standalone!(net::Network, filename::String)
  net.save_standalone = true
  net.standalone_filename = filename
end

end # module
