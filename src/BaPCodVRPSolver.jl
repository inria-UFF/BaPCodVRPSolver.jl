module BaPCodVRPSolver

using JuMP, LightGraphs, Printf
include("blockdecomposition/BlockDecomposition.jl")
include("blockdecompositionextras/BlockDecompositionExtras.jl")
include("bapcod/BaPCod.jl")
include("rcsp/RCSP.jl")
using .BaPCod, .BlockDecomposition, .BlockDecompositionExtras, .RCSP

import Base.show

export VrpModel, VrpGraph, VrpOptimizer
export add_resource!,
       set_resource_bounds!,
       add_elem_set_to_vertex_init_ng_neighbourhood!,
       add_arc!,
       add_arc_var_mapping!,
       set_arc_consumption!,
       set_arc_resource_bounds!,
       get_arc_consumption,
       add_elem_set_to_arc_init_ng_neighbourhood!,
       get_arc_set,
       set_arc_packing_sets!,
       set_vertex_packing_sets!,
       set_additional_arc_elementarity_sets!,
       set_additional_vertex_elementarity_sets!,
       add_graph!,
       define_elementarity_sets_distance_matrix!,
       add_capacity_cut_separator!,
       add_strongkpath_cut_separator!,
       set_branching_priority!,
       enable_rank1_cuts!,
       disable_rank1_cuts!,
       enable_resource_consumption_branching!,
       enable_packset_ryanfoster_branching!,
       optimize!,
       set_cutoff!,
       get_objective_value,
       get_value,
       get_values,
       get_number_of_positive_paths,
       add_cut_callback!,
       add_dynamic_constr!,
       show,
       get_complete_formulation,
       print_enum_paths,
       add_permanent_ryanfoster_constraint!

       @deprecate add_resource(graph; main = false, binary = false, disposable = true, step_size = 0.0) add_resource!(graph; main = false, binary = false, disposable = true, step_size = 0.0)
       @deprecate set_resource_bounds(graph, vertex, resource_id, lb, ub) set_resource_bounds!(graph, vertex, resource_id, lb, ub)
       @deprecate add_pack_set_to_vertex_init_ng_neighbourhood!(graph, vertex, ps_id) add_elem_set_to_vertex_init_ng_neighbourhood!(model, graph, vertex, es_id)
       @deprecate add_arc(graph, tail, head, vars) add_arc!(graph, tail, head, vars)
       @deprecate add_arc_var_mapping(graph, arc_id, vars) add_arc_var_mapping!(graph, arc_id, vars)
       @deprecate set_arc_consumption!(graph, arc_id, res_id, value) set_arc_consumption!(graph, arc_id, res_id, value)
       @deprecate add_pack_set_to_arc_init_ng_neighbourhood!(graph, arc_id, ps_id) add_elem_set_to_arc_init_ng_neighbourhood!(model, graph, arc_id, es_id)
       @deprecate set_arc_packing_sets(model, collection) set_arc_packing_sets!(model, collection)
       @deprecate set_vertex_packing_sets(model, collection) set_vertex_packing_sets!(model, collection)
       @deprecate add_graph(model, graph) add_graph!(model, graph)
       @deprecate define_packing_sets_distance_matrix!(model, matrix) define_elementarity_sets_distance_matrix!(model, graph, matrix)
       @deprecate add_capacity_cut_separator(model, demands, capacity) add_capacity_cut_separator!(model, demands, capacity)
       @deprecate set_branching_priority(model, target, priority) set_branching_priority!(model, target, priority)
       @deprecate enable_rank1_cuts(model) enable_rank1_cuts!(model)
       @deprecate disable_rank1_cuts(model) disable_rank1_cuts!(model)
       @deprecate enable_resource_consumption_branching(model, priority) enable_resource_consumption_branching!(model, priority)
       @deprecate optimize(optimizer) optimize!(optimizer)
       @deprecate set_cutoff(optimizer, value) set_cutoff!(optimizer, value)
       @deprecate add_cut_callback(model, callback, name) add_cut_callback!(model, callback, name)
       @deprecate add_dynamic_constr(optimizer, vars, coeffs, sense, rhs, name) add_dynamic_constr!(optimizer, vars, coeffs, sense, rhs, name)

@enum SetType NoSet=0 ArcSet=1 VertexSet=2

mutable struct VrpVertex
   user_id::Int
   id::Int
   packing_set::Int
   elem_set::Int
   res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for binary resources
   ng_set::Array{Int}
end

mutable struct VrpArc
   id::Int
   tail::Int
   head::Int
   packing_set::Int
   elem_set::Int
   res_consumption::Array{Float64}
   res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for non-binary resources
   vars::Array{Tuple{JuMP.Variable, Float64}}
   ng_set::Array{Int}
end

mutable struct VrpResource
   id::Int
   is_main::Bool
   is_binary::Bool
   is_disposable::Bool
   is_automatic::Bool
   step_size::Float64
end

mutable struct VrpGraph
   id::Int
   source_id::Int
   sink_id::Int
   vertices::Vector{VrpVertex}
   arcs::Vector{VrpArc}
   arc_id_to_arc::Dict{Int,VrpArc}
   incoming_arcs::Vector{Vector{VrpArc}}
   resources::Array{VrpResource}
   multiplicity::Tuple{Int,Int}
   user_vertex_id_map::Dict{Int,Int}
   cycle_problem::Bool # true when source=sink for the user (vertices[sink_id] is only an internal vertex)
   res_bounds_vertex::Dict{Tuple{Int,Int},Tuple{Float64,Float64}} # stores all intervals defined by the user on vertices
   res_bounds_arc::Dict{Tuple{Int,Int},Tuple{Float64,Float64}} # stores all intervals defined by the user on arcs
   arc_rcsp_id_to_id::Dict{Int,Int}
   net::Union{Network,Nothing}
   es_dist_matrix
   elem_sets::Array{Array{Int,1}, 1}
end

mutable struct CapacityCutInfo
   demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}
   capacity::Float64
end

mutable struct VrpModel
   formulation::JuMP.Model
   graphs::Array{VrpGraph}
   packing_sets::Array{Array{Tuple{VrpGraph,Int}, 1}, 1}
   packing_sets_type::SetType
   elem_sets_type::SetType
   define_covering_sets::Bool
   branching_priorities::Dict{String,Int}
   branching_exp_families::Array{Any}
   branching_exps::Array{Any}
   use_rank1_cuts::Bool
   optimizer::Any
   callbacks::Dict{String,Any}
   cap_cuts_info::Array{CapacityCutInfo}
   strongkpath_cuts_info::Array{CapacityCutInfo}
   arcs_by_packing_set_pairs::Array{Array{Tuple{VrpGraph, VrpArc}}, 2}
   ryanfoster_constraints::Vector{Tuple{Integer,Integer,Bool}} # (firstPackSetId,secondPackSetId,together)
   save_standalone::String
end

mutable struct DynamicConstrInfo
   vars::Array{Any}
   coeffs::Array{Float64}
   sense::Any
   rhs::Float64
end

mutable struct CallbackInfo
   constr_name::String
   sep_function::Any
   aux_function::Any
   to_add_constrs::Array{DynamicConstrInfo}
   nb_added_constrs::Int
   cb_idx::Int
end

mutable struct VrpOptimizer
   user_model::VrpModel
   param_file::String
   instance_name::String
   formulation::JuMP.Model
   user_var_to_vars::Dict{JuMP.Variable,Array{JuMP.Variable}}
   callbacks::Dict{String,CallbackInfo}
   cols_in_sol::Array{Tuple{Float64, Dict{JuMP.Variable,Float64}}}
   integer_objective::Bool
   initUB::Float64
   mapped_container_names::Vector{String}
   ignored_container_names::Vector{String}
   stats::Dict{} # execution statistics
   baptreedot_file::String
end

contains(p, s) = findnext(s, p, 1) != nothing

function check_id(id::Int, min_id::Int, max_id::Int)
   (id < min_id || id > max_id) && error("Unknown id $id.")
end

function check_vertex_id(graph::VrpGraph, id::Int)
   !(id in keys(graph.user_vertex_id_map)) && error("Unknown vertex id $id.")
end

function resource_id_in_bapcod(res::VrpResource, graph::VrpGraph)
   if !res.is_binary
      return res.id
   end
   b_id = 1
   for res2 in graph.resources
      if res2.id == res.id
         return b_id
      elseif res2.is_binary
         b_id += 1
      end
   end
end

"""
    VrpModel()

Create an empty VRPSolver model.

It is the main object of the VRPSolver. It is responsible to keep the RCSP graphs (of type VrpGraph), formulation, packing sets, rounded capacity cuts and other definitions (like branching).

"""
function VrpModel(;save_standalone="")
   VrpModel(VrpBlockModel(), VrpGraph[], 
           Array{Tuple{VrpGraph,Int}, 1}[], NoSet, NoSet, false,
           Dict{String,Int}(), Any[], Any[], true,
	   nothing, Dict{String,Any}(), CapacityCutInfo[], CapacityCutInfo[], Array{Array{Tuple{VrpGraph, VrpArc}},2}(undef, 0, 0),
      Vector{Tuple{Integer,Integer,Bool}}(), save_standalone)
end

"""
    enable_rank1_cuts!(model::VrpModel)

Enable use of *limited memory rank-1* cuts during the execution.

These cuts depends of the collection of packing sets (whether in arcs or vertices).
Therefore, it will not be applied if you do not define packing sets.
These cuts are potentially very strong, but each added cut can make the pricing subproblems harder.

"""
function enable_rank1_cuts!(model::VrpModel)
   model.use_rank1_cuts = true
end

function disable_rank1_cuts!(model::VrpModel)
   model.use_rank1_cuts = false
end

"""
    VrpGraph(model::VrpModel, nodes::Array{Int}, source::Int, sink::Int, multiplicity::Tuple{Int,Int})

Create a graph for pricing.

In the VRPSolver, the pricing subproblem is modeled as a *Resource Constrained Shortest Problem* (RCSP). 
VRPSolver will produce paths for the RCSP over the `VrpGraph` to generate the columns.  
Each `VrpGraph` has special `source` and `sink` vertices. The source and sink may be distinct vertices, but may also be the same vertex.

This function does not create a complete `VrpGraph`, it is necessary to create the arcs, a set of resources and to define the resource consumption intervals. 
In addition, it is necessary to map some formulation variables to the arcs of this graph.

See also: [`add_arc!`](@ref), [`add_resource!`](@ref), [`set_resource_bounds!`](@ref), [`set_arc_resource_bounds!`](@ref), [`set_arc_consumption!`](@ref), [`add_arc_var_mapping!`](@ref)

# Arguments
- `nodes::Array{Int}`: list of nodes (id's) of the `VrpGraph` 
- `source::Int`: id of the source node
- `sink::Int`: id of the sink node
- `multiplicity::Tuple{Int,Int}`: multiplicity of the pricing subproblem, i.e., is given lower and upper bounds, respectively, on the number of paths from this graph in a solution.  

# Examples
```julia
# let `model` be a VrpModel
# Creates a VrpGraph with 5 nodes, where at most 2 paths in this graph may appear in the solution
graph1 = VrpGraph(model, [1,2,3,4,5], 1, 5, (0,2)) # paths must start at 1 and end at 5
graph2 = VrpGraph(model, [1,2,3,4,5], 1, 1, (0,2)) # another graph whose paths must start and end at 1 (cycles)
```
"""
function VrpGraph(model::VrpModel, nodes::Array{Int}, source::Int, sink::Int, multiplicity::Tuple{Int,Int})

   vertices = [VrpVertex(nodes[i], i, -1, -1, Dict{Float64, Float64}(), Int[]) for i in 1:length(nodes)]
   user_vertex_id_map = Dict{Int,Int}()
   for i in 1:length(nodes)
      user_vertex_id_map[nodes[i]] = i
   end

   network_source, network_sink = user_vertex_id_map[source], user_vertex_id_map[sink]

   cycle_problem = false
   if source == sink # create a sink node internally?
      network_sink = length(nodes)+1
      push!(vertices, VrpVertex(maximum(nodes)+1, network_sink, -1, -1, Dict{Int,Tuple{Float64, Float64}}(), Int[]))
      cycle_problem = true
   end

   incoming_arcs = [VrpArc[] for i in 1:length(vertices)]
   VrpGraph(-1, network_source, network_sink, vertices, VrpArc[], Dict{Int,VrpArc}(), incoming_arcs, VrpResource[],
         multiplicity, user_vertex_id_map, cycle_problem, Dict{Tuple{Int,Int}, Tuple{Float64, Float64}}(), 
         Dict{Tuple{Int,Int}, Tuple{Float64, Float64}}(), Dict{Int,Int}(), nothing, nothing,
         Array{Array{Int,1}, 1}[]
      )
end

"""
    add_graph!(model::VrpModel, graph::VrpGraph)

Add `VrpGraph` to a `VrpModel`.

This function should be called once for each graph in the model.
"""
function add_graph!(model::VrpModel, graph::VrpGraph)
   graph.id = length(model.graphs) + 1
   push!(model.graphs, graph)

   #first, create automatic resource if needed
   if isempty(filter(r -> r.is_main, graph.resources))
      res_id = add_resource!(graph, main = true, step_size = 1.0)
      graph.resources[res_id].is_automatic = true
      println("Warning: An automatic main resource for Graph $(graph.id) has been created")
      flush(stdout)
   end

   # second, define intervals on vertices
   for key in keys(graph.res_bounds_vertex)
      vertex, res_id = key
      res = graph.resources[res_id] 
      if res.is_binary
         if vertex != graph.source_id
            graph.vertices[vertex].res_bounds[res_id] = graph.res_bounds_vertex[key]
         end
      else
         for arc in graph.incoming_arcs[vertex]
            arc.res_bounds[res_id] = graph.res_bounds_vertex[key]
         end
      end
   end

   # third, define intervals on arcs (overriding, if necessary, the previous definitions)
   for key in keys(graph.res_bounds_arc)
      arc_id, res_id = key
      graph.arcs[arc_id].res_bounds[res_id] = graph.res_bounds_arc[key]
   end

   return graph
end

"""
    add_resource!(graph::VrpGraph; main = false, binary = false, disposable = true, step_size = 0.0)

Add a resource to the VrpGraph `graph`. 

# Optional arguments 
- `main::Bool`: indicates that the resource is main or secondary. 
- `binary::Bool`: indicates that the resource is binary, i.e., its accumulated consumption can only be `0` or `1`. 
- `disposable::Bool`: indicates that the resource is disposable or non-disposable.
- `step_size::Float64`: only for main resources, this advanced parameter is used for determining the resource consumption intervals that define each bucket on the labeling algorithm during the pricing. 
if step_size is not given, it is determined automatically for main resources, based on the parameter `RCSPnumberOfBucketsPerVertex`.

# Examples
```julia
# let `graph` be a VrpGraph
r1 = add_resource!(graph) # create a secondary, disposable, non-binary resource
r2 = add_resource!(graph, main=true, disposable=false) # create a main, non-disposable, non-binary resource
```

"""
function add_resource!(graph::VrpGraph; main = false, binary = false, disposable = true, step_size = 0.0)
   main && binary && error("VRPSolver error: binary resource cannot be main resource")

   res_id = length(graph.resources) + 1
   push!(graph.resources, VrpResource(res_id, main, binary, disposable, false, step_size))
   if binary
      for vertex in graph.vertices
         vertex.res_bounds[res_id] = (0.0, 1.0)
      end
   end
   for arc in graph.arcs
      push!(arc.res_consumption, 0.0)
      if !binary
         arc.res_bounds[res_id] = (0.0, 0.0)
      end
   end
   return res_id
end

function set_resource_bounds_aux!(graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64)
   if !graph.cycle_problem && vertex == graph.source_id && (lb != 0.0 || ub != 0.0)
      println("VRPSolver warning: Interval set for resource $(res_id) on source node (when source != sink) is ignored (by definition is [0.0,0.0])")
      flush(stdout)
   end
   if graph.resources[res_id].is_binary
      if (lb,ub) != (0,0) && (lb,ub) != (0,1) && (lb,ub) != (1,1)
         error("VRPSolver error: binary resources only supports the intervals [0,0], [0,1] and [1,1].")
      end
      if !graph.resources[res_id].is_disposable && ((vertex == graph.sink_id) || (graph.cycle_problem && (vertex == graph.source_id))) && (lb != ub)
         error("VRPSolver error: non-disposable binary resources cannot have the interval [0,1] set for the sink node")
      end
   end
   graph.res_bounds_vertex[vertex,res_id] = (lb, ub) # store interval for res_id on vertex (to be defined in add_graph)
   if graph.cycle_problem && (vertex == graph.source_id) # add bounds for graph.sink_id (because the user will not set)
      set_resource_bounds_aux!(graph, graph.sink_id, res_id, lb, ub)
   end
end

"""
    set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64)

Set the resource bounds to a vertex of the VrpGraph `graph`.

Defining the interval ``[lb,ub]`` for res_id at vertex is equivalent to defining the same interval for every incoming arc of vertex.

# Arguments
- `vertex::Int`: vertex id in the VrpGraph `graph`.
- `res_id::Int`: resource id in the VrpGraph `graph`.
- `lb::Float64`: lower bound on the resource consumption at the vertex.
- `ub::Float64`: upper bound on the resource consumption at the vertex.

"""
function set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64)
   check_id(res_id, 1, length(graph.resources))
   set_resource_bounds_aux!(graph, graph.user_vertex_id_map[vertex], res_id, lb, ub)
end

set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb, ub) = set_resource_bounds!(graph, vertex, res_id, Float64(lb), Float64(ub))

"""
    add_elem_set_to_vertex_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int)

Add elementarity set to the vertex ng-set. 

This is an explicit way to set initial ng-neighbourhoods, which has priority over the definition with [`define_elementarity_sets_distance_matrix!`](@ref).
In fact, distance matrix is taken into account only if user-defined ng-sets are empty.

# Arguments
- `model::VrpModel`: model
- `graph::VrpGraph`: graph
- `vertex_id::Int`: vertex id
- `es_id::Int`: elementarity set id. The valid ids are ``\\{1,2,\\dots,|\\mathcal{P}^V|,|\\mathcal{P}^V|+1,|\\mathcal{P}^V|+2,\\dots,|\\mathcal{P}^V|+|\\mathcal{E}^V|\\}`` (considering automatic and additional elem. sets)
"""
function add_elem_set_to_vertex_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int)
   check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
   check_vertex_id(graph, vertex_id)
   vertex = graph.vertices[graph.user_vertex_id_map[vertex_id]]
   if !(es_id in vertex.ng_set)
      push!(vertex.ng_set, es_id)
   end
end

"""
    add_elem_set_to_arc_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int)

Add elementarity set to the arc ng-set.

This is an explicit way to set initial ng-neighbourhoods, which has priority over the definition with [`define_elementarity_sets_distance_matrix!`](@ref).
In fact, distance matrix is taken into account only if user-defined ng-sets are empty.


# Arguments
- `model::VrpModel`: model
- `graph::VrpGraph`: graph
- `arc_id::Int`: arc id
- `es_id::Int`: elementarity set id. The valid ids are ``\\{1,2,\\dots,|\\mathcal{P}|,|\\mathcal{P}|+1,|\\mathcal{P}|+2,\\dots,|\\mathcal{P}|+|\\mathcal{E}|\\}`` (considering automatic and additional elem. sets)
"""
function add_elem_set_to_arc_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int)
   check_id(arc_id, 1, length(graph.arcs))
   check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
   arc = graph.arcs[arc_id]
   if !(es_id in arc.ng_set)
      push!(arc.ng_set, es_id)
   end
end

#TODO: check for repeated variables before passing the model to bapcod
"""
    add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, vars::Array{Tuple{JuMP.Variable, Float64}})

Define variable mapping for an existing arc.

# Arguments

- `arc_id::Int`: id of the arc
- `vars::Array{Tuple{JuMP.Variable, Float64}}`: variables to be mapped to the arc. It is a set of pairs of variable and coefficient. There are extensions where `vars` is `Array{JuMP.Variable}` and `JuMP.Variable` which consider all coefficients as `1.0`.

# Examples
```julia
add_arc_var_mapping!(graph, arc_id, [(x1,2.0), (x2,2.0)]) # map 2x1 and 2x2 to the arc in `graph` with id 1
add_arc_var_mapping!(graph, arc_id, [x1, x2]) # map x1 and x2
add_arc_var_mapping!(graph, arc_id, x) # map x
```
"""
function add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, vars::Array{Tuple{JuMP.Variable, Float64}})
   check_id(arc_id, 1, length(graph.arcs))
   for (user_var, coeff) in vars
      if getcategory(user_var) == :Bin
         error("ERROR: mapping a binary variable is not allowed. Please redefine it as integer or continuous.")
      end
      for (var, c) in graph.arcs[arc_id].vars
         if var == user_var
	    error("ERROR: variable $(user_var) is mapped more than once to arc $(arc_id) of graph $(graph.id)")
	 end
      end
      push!(graph.arcs[arc_id].vars, (user_var, coeff))
   end
end
add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, var::JuMP.Variable) = add_arc_var_mapping!(graph, arc_id, [(var, 1.0)])
add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, vars::Array{JuMP.Variable}) = add_arc_var_mapping!(graph, arc_id, [(var, 1.0) for var in vars])

default_consumptions(graph::VrpGraph) = [0.0 for i in 1:length(graph.resources)]  

function default_resbounds(graph::VrpGraph)
   bounds = Dict{Int,Tuple{Float64,Float64}}()
   for res in graph.resources  
      if !res.is_binary
         bounds[res.id] = (0.0, 0.0)
      end
   end
   return bounds
end

"""
    add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{Tuple{JuMP.Variable, Float64}} = Tuple{JuMP.Variable, Float64}[])

Add arc `(tail,head)` to `graph` and return the arc id. 

Adding parallel arcs is allowed, since they will have different identifiers in `graph`.

# Optional argument
- `vars::Array{Tuple{JuMP.Variable, Float64}}`: variables to be mapped to the arc `(tail,head)`. It is a set of pairs of variable and coefficient. There are extensions where `vars` is `Array{JuMP.Variable}` and `JuMP.Variable` which consider all coefficients as `1.0`. 

# Examples
```julia
# let `x1` and `x2` two decision variables
add_arc!(graph, 1, 2, [(x1,1.0), (x2,2.0)]) # add (1,2) mapped to x1 with coefficient 1 and to x2 with coefficient 2 
add_arc!(graph, 2, 1, [x1, x2]) # add (2,1) mapped to x1 and x2
arc_id = add_arc!(graph, 1, 2, x1) # add (1,2) mapped to x and get the arc id
arc_id = add_arc!(graph, 1, 2) # add (1,2) without mapped variable
```
"""
function add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{Tuple{JuMP.Variable, Float64}} = Tuple{JuMP.Variable, Float64}[])
   arc = []
   if graph.cycle_problem && (graph.user_vertex_id_map[head] == graph.source_id)
      arc = VrpArc(length(graph.arcs) + 1, graph.user_vertex_id_map[tail], graph.sink_id, -1, -1, default_consumptions(graph), default_resbounds(graph), vars, Int[])
   else
      arc = VrpArc(length(graph.arcs) + 1, graph.user_vertex_id_map[tail], graph.user_vertex_id_map[head], -1, -1, default_consumptions(graph), default_resbounds(graph), vars, Int[])
   end
   push!(graph.incoming_arcs[arc.head], arc)
   push!(graph.arcs, arc)
   graph.arc_id_to_arc[arc.id] = arc
   return arc.id
end
add_arc!(graph::VrpGraph, tail::Int, head::Int, var::JuMP.Variable) = add_arc!(graph, tail, head, [(var, 1.0)])
add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{JuMP.Variable}) = add_arc!(graph, tail, head, [(var, 1.0) for var in vars])

"""
    set_arc_consumption!(graph::VrpGraph, arc_id::Int, res_id::Int, value::Union{Int,Float64})

Set the arc consumption for a specific resource.

# Arguments
- `graph::VrpGraph`: graph to be considered
- `arc_id::Int`: arc to be considered
- `res_id::Int`: resource id to define consumption 
- `value::Union{Int,Float64}`: consumption value which can be integer or real

# Example
```julia
set_arc_consumption!(graph, 3, 1, 2.5) # define a consumption of 2.5 for the resource 1 when passing by the arc 3 
```
"""
function set_arc_consumption!(graph::VrpGraph, arc_id::Int, res_id::Int, value::Union{Int,Float64})
   check_id(arc_id, 1, length(graph.arcs))
   graph.resources[res_id].is_binary && (value != -1 && value != 0 && value != 1) && error("VRPSolver error: arc consumption for binary resources must be -1, 0, or 1")
   graph.arcs[arc_id].res_consumption[res_id] = Float64(value)
end

"""
    get_arc_consumption(graph::VrpGraph, arc_id::Int, res_id::Int)

Get the arc consumption value for a specific resource.
"""
function get_arc_consumption(graph::VrpGraph, arc_id::Int, res_id::Int)
   graph.arcs[arc_id].res_consumption[res_id]
end

"""
    set_arc_resource_bounds!(graph::VrpGraph, arc_id::Int, res_id::Int, lb::Float64, ub::Float64)

Set the resource bounds for an arc of the VrpGraph `graph`.

# Arguments
- `arc::Int`: arc id in the VrpGraph `graph`.
- `res_id::Int`: resource id in the VrpGraph `graph`.
- `lb::Float64`: lower bound on the resource consumption at the vertex.
- `ub::Float64`: upper bound on the resource consumption at the vertex.

"""
function set_arc_resource_bounds!(graph::VrpGraph, arc_id::Int, res_id::Int, lb::Float64, ub::Float64)
   check_id(res_id, 1, length(graph.resources))
   res = graph.resources[res_id] 
   if res.is_binary
      error("ERROR: Resource bounds on arcs for binary resources is not yet implemented. Please use the vertex-based function set_resource_bounds!.")
   else
      graph.res_bounds_arc[arc_id,res_id] = (lb, ub)
   end
end

set_arc_resource_bounds!(graph::VrpGraph, arc_id::Int, res_id::Int, lb, ub) = set_arc_resource_bounds!(graph, arc_id, res_id, Float64(lb), Float64(ub))

function is_preprocessed_arc(graph::VrpGraph, arc::VrpArc)
   if !graph.cycle_problem &&
        (arc.head == graph.source_id || arc.tail == graph.sink_id)
     return true
   end
   return false
end

"""
    get_arc_set(graph::VrpGraph)

Return the set of arcs of a VrpGraph as an array of triples, where each triple contains the arc (first two elements) and its id (third element).

For example, if the function returns `[(1,2,1), (2,1,2)]`, it must be interpreted as `graph` having two arcs: `(1,2)` with id `1` and `(2,1)` with id `2`.
"""
function get_arc_set(graph::VrpGraph)
   arcs = []
   if graph.cycle_problem
      for arc in graph.arcs
         i = (arc.tail == graph.sink_id) ? graph.source_id : arc.tail
         j = (arc.head == graph.sink_id) ? graph.source_id : arc.head
         push!(arcs, (graph.vertices[i].user_id, graph.vertices[j].user_id, arc.id))
      end
   else
      arcs = [(graph.vertices[arc.tail].user_id, graph.vertices[arc.head].user_id, arc.id) for arc in graph.arcs]
   end
   return arcs
end

function reset_packing_sets(user_model::VrpModel)
   empty!(user_model.packing_sets)
   user_model.packing_sets_type = NoSet
   for graph in user_model.graphs
      for vertex in graph.vertices
         vertex.packing_set = -1
         empty!(vertex.ng_set)
      end
      for arc in graph.arcs
         arc.packing_set = -1
         empty!(arc.ng_set)
      end
   end
end

function reset_elem_sets(user_model::VrpModel)
   user_model.elem_sets_type = NoSet
   for graph in user_model.graphs
      empty!(graph.elem_sets)
      graph.es_dist_matrix = nothing
      for vertex in graph.vertices
         vertex.elem_set = -1
         empty!(vertex.ng_set)
      end
      for arc in graph.arcs
         arc.elem_set = -1
         empty!(arc.ng_set)
      end
   end
end

function add_arc_to_packing_set(model::VrpModel, graph::VrpGraph, arc_id::Int, packing_set_id::Int)
   check_id(arc_id, 1, length(graph.arcs))
   check_id(packing_set_id, 1, length(model.packing_sets))
   (graph.arcs[arc_id].packing_set != -1) && error("an arc cannot belong to more than 1 packing set")
   graph.arcs[arc_id].packing_set = packing_set_id
   (graph.arcs[arc_id].elem_set != -1) && error("an arc cannot belong to more than 1 elementarity set")
end

function add_arc_to_elementarity_set(model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int)
   check_id(arc_id, 1, length(graph.arcs))
   check_id(es_id, 1, length(graph.elem_sets) + length(model.packing_sets))
   (graph.arcs[arc_id].elem_set != -1) && error("an arc cannot belong to more than 1 elementarity set")
   (graph.arcs[arc_id].packing_set != -1) && error("an arc cannot belong to more than 1 elementarity set")
   graph.arcs[arc_id].elem_set = es_id
end

"""
    set_arc_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})

Define a collection of packing sets on arcs. For each defined packing set, VRPSolver automatically will create an equivalent elementarity set on arcs.

The collection must be a set of mutually disjoint subsets of arcs. Not all arcs need to belong to some packing set.
The index of a packing set in the array defines its packing set id.

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# A packing set is a list of pairs of VrpGraph and arc id
ps_1 = [(G[1],1),(G[2],1),(G[3],1)] # packing set composed by arcs of different graphs
ps_2 = [(G[2],2),(G[2],4)] # packing set composed by arcs of the same graph
ps_3 = [(G[2],3),(G[3],2)] # another packing set
set_arc_packing_sets!(model, [ps_1, ps_2, ps_3]) # passing the collection of packing sets to the model
```
"""
function set_arc_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})
   if user_model.elem_sets_type != NoSet
      error("Packing sets cannot be defined after elementarity sets")
   end
   reset_packing_sets(user_model)
   user_model.packing_sets = collection
   user_model.packing_sets_type = ArcSet
   for ps_id in 1:length(collection)
      for (graph, arc_id) in collection[ps_id]
         add_arc_to_packing_set(user_model, graph, arc_id, ps_id)
      end
   end
end

"""
    set_additional_arc_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})

Define an additional collection of elementarity sets on arcs.

The collection must be a set of mutually disjoint subsets of arcs. Not all arcs need to belong to some elementarity set.
The index ``i`` of a elementarity set in the array defines its elementarity set id as ``|\\mathcal{P}|+i`` because there are ``|\\mathcal{P}|`` 
automatic elementarity sets (one for each packing set) created with [`set_arc_packing_sets!`](@ref).

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# An elementarity set is a pair of VrpGraph and a list of arc ids
es_1 = (G[1],[1,2,4]) # elem. set on graph G[1] composed by the arcs 1, 2, and 4  
es_2 = (G[2],[2,3,5,9]) # elem. set on G[2] with 4 arcs
es_3 = (G[2],[1,4,6,7]) # another elem. set on G[2]
set_additional_arc_elementarity_sets!(model, [es_1, es_2, es_3]) # passing the collection of elem. sets to the model
```
"""
function set_additional_arc_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})
   if user_model.packing_sets_type == VertexSet
      error("Vertex packing sets and arc elementarity sets are not compatible")
   end
   reset_elem_sets(user_model)
   user_model.elem_sets_type = ArcSet
   for es_id in 1:length(collection)
      graph, arc_ids = collection[es_id]
      push!(graph.elem_sets, arc_ids)
      for arc_id in arc_ids
         add_arc_to_elementarity_set(user_model, graph, arc_id, es_id)
      end
   end
   return length(user_model.packing_sets)
end

function add_vertex_to_packing_set(model::VrpModel, graph::VrpGraph, vertex_id::Int, packing_set_id::Int)
   check_vertex_id(graph, vertex_id)
   check_id(packing_set_id, 1, length(model.packing_sets))
   vertexAlgId = graph.user_vertex_id_map[vertex_id]
   (graph.vertices[vertexAlgId].packing_set != -1) && error("a vertex cannot belong to more than 1 packing set")
   graph.vertices[vertexAlgId].packing_set = packing_set_id
   (graph.vertices[vertexAlgId].elem_set != -1) && error("a vertex cannot belong to more than 1 elementarity set")
end

function add_vertex_to_elem_set(model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int)
   check_vertex_id(graph, vertex_id)
   check_id(es_id, 1, length(graph.elem_sets) + length(model.packing_sets))
   vertexAlgId = graph.user_vertex_id_map[vertex_id]
   (graph.vertices[vertexAlgId].elem_set != -1) && error("a vertex cannot belong to more than 1 elementarity set")
   graph.vertices[vertexAlgId].elem_set = es_id
   (graph.vertices[vertexAlgId].packing_set != -1) && error("a vertex cannot belong to more than 1 elementarity set")
end

"""
    set_vertex_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})

Define a collection of packing sets on vertices. For each defined packing set, VRPSolver automatically will create an equivalent elementarity set on vertices.

The collection must be a set of mutually disjoint subsets of vertices. Not all vertices need to belong to some packing set.
The index of a packing set in the array defines its packing set id.

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# A packing set is a list of pairs of VrpGraph and vertex id
ps_1 = [(G[1],1),(G[2],1),(G[3],1)] # packing set composed by vertices of different graphs
ps_2 = [(G[2],2),(G[2],4)] # packing set composed by vertices of the same graph
ps_3 = [(G[2],3),(G[3],2)] # another packing set
set_vertex_packing_sets!(model, [ps_1, ps_2, ps_3]) # passing the collection of packing sets to the model
```
"""
function set_vertex_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1}, define_covering_sets::Bool=false)
   if user_model.elem_sets_type != NoSet
      error("Packing sets cannot be defined after elementarity sets")
   end
   reset_packing_sets(user_model)
   user_model.packing_sets = collection
   user_model.packing_sets_type = VertexSet
   user_model.define_covering_sets = define_covering_sets
   n = length(collection)
   for ps_id in 1:n
      for (graph, vertex_id) in collection[ps_id]
         add_vertex_to_packing_set(user_model, graph, vertex_id, ps_id)
      end
   end

   # function to compute the packing set pair connected by the arc of a pair
   # (graph, arc) where vectices not associated to packing sets are assigned to
   # the dummy packing set id n+1
   function get_packing_set_pair(gr_arc::Tuple{VrpGraph, VrpArc})::Tuple{Int, Int}
      graph = gr_arc[1]
      arc = gr_arc[2]
      head = graph.vertices[arc.head].packing_set
      head = (head < 1) ? n + 1 : head
      tail = graph.vertices[arc.tail].packing_set
      tail = (tail < 1) ? n + 1 : tail
      if head < tail
         return head, tail
      else
         return tail, head
      end
   end

   # build data structures to access lists of arcs by packing set pairs
   # and by mapped variables (to be used next)
   user_model.arcs_by_packing_set_pairs = [
      Tuple{VrpGraph, VrpArc}[] for i in 1:n, j in 1:n
   ]
   mapped_arcs_by_vars = Dict{JuMP.Variable, Array{Tuple{VrpGraph, VrpArc}, 1}}()
   for graph in user_model.graphs
      for arc in graph.arcs
         head, tail = get_packing_set_pair((graph, arc))
         if tail <= n            
            push!(user_model.arcs_by_packing_set_pairs[head, tail], (graph, arc))
         end
         for (var, val) in arc.vars
            if haskey(mapped_arcs_by_vars, var)
               push!(mapped_arcs_by_vars[var], (graph, arc))
            else
               mapped_arcs_by_vars[var] = [(graph, arc)]
            end
         end
      end
   end

   # determine, for each packing set pairs, which arcs are not covered by
   # appropriate variables (variables mapped to arcs connecting the same packing
   # set pair, in any direction)
   for (var, gr_arcs) in mapped_arcs_by_vars
      head, tail = get_packing_set_pair(gr_arcs[1])
      appropriate = (head != tail) && (tail <= n)
      for gr_arc in gr_arcs[2:end]
         if (head, tail) != get_packing_set_pair(gr_arc)
            appropriate = false
         end
      end
      if appropriate
         ps_gr_arcs = user_model.arcs_by_packing_set_pairs[head, tail]
         filter!(x -> !(x in gr_arcs), ps_gr_arcs)
      end
   end
end

"""
    set_additional_vertex_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})

Define an additional collection of elementarity sets on vertices.

The collection must be a set of mutually disjoint subsets of vertices. Not all vertices need to belong to some elementarity set.
The index ``i`` of a elementarity set in the array defines its elementarity set id as ``|\\mathcal{P}^V|+i`` because there are ``|\\mathcal{P}^V|`` 
automatic elementarity sets (one for each packing set) created with [`set_vertex_packing_sets!`](@ref).

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# An elementarity set is a pair of VrpGraph and a list of vertex ids
es_1 = (G[1],[1,2,4]) # elem. set on graph G[1] composed by the vertices 1, 2, and 4  
es_2 = (G[2],[2,3,5,9]) # elem. set on G[2] with 4 vertices
es_3 = (G[2],[1,4,6,7]) # another elem. set on G[2]
set_additional_vertex_elementarity_sets!(model, [es_1, es_2, es_3]) # passing the collection of elem. sets to the model
```
"""
function set_additional_vertex_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})
   if user_model.packing_sets_type == ArcSet
      error("Arc packing sets and vertex elementarity sets are not compatible")
   end
   reset_elem_sets(user_model)
   user_model.elem_sets_type = VertexSet
   n = length(collection)
   for es_id in 1:n
      (graph, vertex_ids) = collection[es_id]
      push!(graph.elem_sets, vertex_ids)
      for vertex_id in vertex_ids
         add_vertex_to_elem_set(user_model, graph, vertex_id, es_id)
      end 
   end
end

"""
    set_branching_priority!(model::VrpModel, var_container_name::String, priority::Int)

Set the branching priority for a decision variable. 

Branching is performed on a variable with priority ``k`` only when there is no possible branching for variables with priority higher than ``k``.

# Arguments
- `var_container_name::String`: JuMP decision variable name.
- `priority::Int`: the priority of the variable, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
@variable(model.formulation, x[i in 1:10] >= 0, Int)
@variable(model.formulation, y[i in 1:20] >= 0, Int)
... 
set_branching_priority!(model, "x", 2) # x has higher priority
set_branching_priority!(model, "y", 1) # than y
```
"""
function set_branching_priority!(model::VrpModel, var_container_name::String, priority::Int)
   model.branching_priorities[var_container_name] = priority
end

"""
    function set_branching_priority!(model::VrpModel, exps_family, name::String, priority::Int)

Set the branching priority for a set of JuMP expressions.

# Arguments
- `exps_family`: set of expressions created with JuMP macro @expression 
- `name::String`: JuMP expression name.
- `priority::Int`: the priority of the expressions, which is higher when this value is increased.

"""
function set_branching_priority!(model::VrpModel, exps_family, name::String, priority::Int)
   push!(model.branching_exp_families, (exps_family, name, priority))
end

"""
    function set_branching_priority!(model::VrpModel, exp::JuMP.GenericAffExpr{Float64,JuMP.Variable}, name::String, priority::Int)

Set the branching priority for a JuMP expression.

# Arguments
- `exp`: expression created with JuMP macro @expression 
- `name::String`: JuMP expression name.
- `priority::Int`: the priority of the expression, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
@variable(model.formulation, x[i in 1:10] >= 0, Int)
@expression(model.formulation, exp, sum(x[i] for i in 1:5))
... 
set_branching_priority!(model, "x", 2) # x has higher priority
set_branching_priority!(model, exp, "exp", 1) # than exp
```
"""
function set_branching_priority!(model::VrpModel, exp::JuMP.GenericAffExpr{Float64,JuMP.Variable}, name::String, priority::Int)
   push!(model.branching_exps, (exp, name, priority))
end

"""
    function enable_resource_consumption_branching!(model::VrpModel, priority::Int)

Enable the accumulated resource consumption branching. 
Given a packing set ``S \\in \\mathcal{P}``, a main resource ``r \\in R^k_M``, and a certain threshold value ``t^*``: in the left child make ``u_{a,r}=t^*``,
for all ``a \\in S``; in the right child make ``l_{a,r}=t^*``, for all ``a \\in S``.
The packing set, the main resource, and the threshold are automatically chosen by VRPSolver during the execution.
The branching is not likely to be complete, in the sense that some fractional ``\\lambda`` solutions can not be eliminated by it. 
However, it does not increase the pricing difficulty and it may work well in practice, postponing (and even avoiding) the use of a branching that makes pricing harder.

# Arguments
- `priority::Int`: the branching priority, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
... 
enable_resource_consumption_branching!(model, 1) 
```
"""
function enable_resource_consumption_branching!(model::VrpModel, priority::Int)
   model.branching_priorities["_res_cons_branching_"] = priority
end

"""
    function enable_packset_ryanfoster_branching!(model::VrpModel, priority::Int)

Enable the Ryan and Foster branching.
Given two distinct packing sets ``S`` and ``S'`` in ``\\mathcal{P}``, let ``P(S,S') \\subseteq P`` be the subset of the paths that contain arcs in both ``S`` and ``S'``.
The branching is over the value of ``\\sum_{p\\in P(S,S')} \\lambda_p``, either 0 or 1.
The packing sets are automatically chosen by VRPSolver during the execution.
This branching is still to be avoided if possible, because it makes the pricing harder. However, using that scheme leads to more balanced search trees.

# Arguments
- `priority::Int`: the branching priority, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
... 
enable_packset_ryanfoster_branching!(model, 2) 
```
"""
function enable_packset_ryanfoster_branching!(model::VrpModel, priority::Int)
   model.branching_priorities["_ryanfoster_branching_"] = priority
end

"""
    function add_permanent_ryanfoster_constraint!(model::VrpModel, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool)

Set a permanent Ryan and Foster constraint for a pair of packing sets.

# Arguments
- `firstPackSetId::Int`: first packing set id.
- `secondPackSetId::Int`: second packing set id.
- `together::Bool`: say whether they should be together in the same path or should not be together.

# Example
```julia
model = VrpModel()
... 
add_permanent_ryanfoster_constraint!(model, 1, 2, true) 
add_permanent_ryanfoster_constraint!(model, 3, 4, false) 
```
"""
function add_permanent_ryanfoster_constraint!(model::VrpModel, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool)
   push!(model.ryanfoster_constraints, (firstPackSetId, secondPackSetId, together))
end

"""
    define_elementarity_sets_distance_matrix!(model::VrpModel, graph::VrpGraph, matrix::Array{Array{Float64,1},1})

Define distance matrix between the elementarity sets for a specific graph.

This is a convenient way of defining the initial ng-sets for pricing over ``\\mathcal{E}``-ng-paths: the ng-set of each vertex (arc) is defined as being the ng-size elementarity sets 
that are closer to its own elementarity set, according to the given distance matrix. The element `matrix[i][j]` represents the distance between the elementarity set ``i`` to the 
elementarity set ``j``.
An index of an elementarity set is defined during its creation: automatic elementarity sets (one for each packing set) has the indexes
 ``\\{1,2,\\dots,|\\mathcal{P}|\\}``, whereas the additional elementarity sets has the indexes ``\\{|\\mathcal{P}|+1,|\\mathcal{P}|+2,\\dots\\}``
  (these indexes consider ``|\\mathcal{P}^V|`` instead of ``|\\mathcal{P}|`` for packing sets on vertices). 
To define the distance between two elementarity sets, 
it is recommended to consider some metric which involves all elements (vertices or arcs) of both elementarity sets.

"""
function define_elementarity_sets_distance_matrix!(model::VrpModel, graph::VrpGraph, matrix::Array{Array{Float64,1},1})
   n = length(model.packing_sets) + length(graph.elem_sets)
   size(matrix, 1) != n && error("wrong matrix dimension")
   for i in 1:n
      size(matrix[i], 1) != n && error("wrong matrix dimension")
   end
   graph.es_dist_matrix = matrix
end

function extract_user_var_to_graphs(user_model::VrpModel)
   var_to_graphs = Dict{JuMP.Variable, Array{Int}}()
   for graph_id in 1:length(user_model.graphs)
      for arc_id in 1:length(user_model.graphs[graph_id].arcs)
         for (var,coeff) in user_model.graphs[graph_id].arcs[arc_id].vars
            if !haskey(var_to_graphs, var)
               var_to_graphs[var] = [graph_id]
            elseif !(graph_id in var_to_graphs[var])
               push!(var_to_graphs[var], graph_id)
            end
         end
      end
   end
   return var_to_graphs
end

function split_var_name(var)
   var_name = getname(var)
   var_container_name = contains(var_name, "[") ? split(var_name, "[")[1] : var_name
   var_id = contains(var_name, "[") ? split(split(var_name, "[")[2], "]")[1] : "1"
   var_id = replace(var_id, r"\(" => s"")
   var_id = replace(var_id, r"\)"=> s"")
   var_id = tuple([parse(Int64, x) for x in split(var_id, ",") if x != ""]...)
   return (var_container_name, var_id)
end

function extract_user_var_to_vars(user_form::JuMP.Model, formulation::JuMP.Model, user_var_to_graphs, mapped_names)
   var_to_vars = Dict{JuMP.Variable,Array{JuMP.Variable}}()
   user_vars = [JuMP.Variable(user_form, i) for i in 1:user_form.numCols]
   for user_var in user_vars
      (var_container_name, var_id) = split_var_name(user_var)
      var_container = formulation[Symbol(var_container_name)]
      if !haskey(user_var_to_graphs, user_var) #master variable
         if !(var_container_name in mapped_names)
            ids = [0]
         else
            ids = []
         end
      else
         ids = user_var_to_graphs[user_var] #subproblem variable
      end
      vars = []
      for id in ids
         push!(vars, var_container[(id, var_id)])
      end
      var_to_vars[user_var] = vars
   end
   return var_to_vars
end

function generate_pricing_networks(user_model::VrpModel, user_var_to_vars::Dict{JuMP.Variable,Array{JuMP.Variable}})
   networks = []
   for graph in user_model.graphs
      nbPackSets = length(user_model.packing_sets)
      nbElemSets = nbPackSets + length(graph.elem_sets)
      nbCovSets = 0
      if user_model.define_covering_sets
         nbCovSets = nbPackSets
      end
      network = Network(length(graph.vertices), source = graph.source_id, sink = graph.sink_id, elemsets = nbElemSets,
                        packsets = nbPackSets, covsets = nbCovSets)
      if user_model.save_standalone != ""
         save_standalone!(network, user_model.save_standalone)
      end
                  
      graph.net = network
      #resources and vertices
      for resource in graph.resources
         if resource.is_main
            addmainresource!(network, resource.id, stepsize = resource.step_size, disposable = resource.is_disposable)
         elseif resource.is_binary
            addspecialresource!(network, disposable = resource.is_disposable)
         else
            addresource!(network, resource.id,  disposable = resource.is_disposable)
         end
         for vertex in graph.vertices
            if resource.is_binary
               res_seq_id = resource_id_in_bapcod(resource, graph)
               vertexspecialconsumptionbounds!(network, vertex.id, res_seq_id, lb = vertex.res_bounds[resource.id][1], ub = vertex.res_bounds[resource.id][2])
            else
               vertexconsumptionbounds!(network, vertex.id, resource.id, lb = -1e12, ub = 1e12)
            end
         end
      end
      #adding vertices to packing and elem sets
      for vertex in graph.vertices
         if vertex.packing_set != -1
            add_vertex_to_packing_set!(network, vertex.id, vertex.packing_set)
            attach_elem_set_to_vertex!(network, vertex.id, vertex.packing_set)
            if user_model.define_covering_sets
               add_vertex_to_covering_set!(network, vertex.id, vertex.packing_set)
            end
         elseif vertex.elem_set != -1
            attach_elem_set_to_vertex!(network, vertex.id, nbPackSets + vertex.elem_set)
         end
         #defining the neighbourhood of the vertex
         for es_id in vertex.ng_set
            add_vertex_to_mem_of_elem_set!(network, vertex.id, es_id)
         end
      end
      #arcs
      for arc in graph.arcs
         if is_preprocessed_arc(graph, arc)
            continue
         end
         arc_rcsp_id = add_edge!(network, arc.tail, arc.head)
         graph.arc_rcsp_id_to_id[arc_rcsp_id] = arc.id
         for resource in graph.resources
            if resource.is_binary
               res_seq_id = resource_id_in_bapcod(resource, graph)
               edgespecialconsumption!(network, arc_rcsp_id, res_seq_id, value = arc.res_consumption[resource.id])
            else
               edgeconsumption!(network, arc_rcsp_id, resource.id, value = arc.res_consumption[resource.id])
               arcconsumptionbounds!(network, arc_rcsp_id, resource.id, lb = arc.res_bounds[resource.id][1] , ub = arc.res_bounds[resource.id][2])
            end
         end
         #adding arc to packing_set
         if arc.packing_set != -1
            add_edge_to_packing_set!(network, arc_rcsp_id, arc.packing_set)
            attach_elem_set_to_edge!(network, arc_rcsp_id, arc.packing_set)
         elseif arc.elem_set != -1
            attach_elem_set_to_edge!(network, arc_rcsp_id, nbPackSets + arc.elem_set)
         end
         #defining the neighbourhood of the arc
         for es_id in arc.ng_set
            add_arc_to_mem_of_elem_set!(network, arc_rcsp_id, es_id)
         end
         #arc variables
         vars = []
         for (user_var, coeff) in arc.vars
            for var in user_var_to_vars[user_var]
               (var_container_name, var_id) = split_var_name(var)
               if graph.id == var_id[1]
                  push!(vars, (var, coeff))
               end
            end
         end
         set_edge_vars!(network, arc_rcsp_id, vars)
      end

      #packing sets distance matrix
      if graph.es_dist_matrix != nothing
         set_elementarity_sets_distance_matrix!(network, graph.es_dist_matrix)
      end

      #ryan and foster constraints
      for (ps1, ps2, tog) in user_model.ryanfoster_constraints
         add_permanent_ryan_foster_constraint!(network, ps1, ps2, tog)
      end

      push!(networks, network)
   end
   return networks
end

"""
    add_cut_callback!(user_model::VrpModel, callback::Any, constr_name::String)

Add user cut callback for separation during the execution.

Cuts must be added through the function [`add_dynamic_constr!`](@ref).

# Arguments
- `model::VrpModel`: model to be added the cut callback.
- `callback::Any`: function with the separation algorithm 
- `constr_name::String`: a nome for the cut callback 

# Example
```julia
# let model be a VrpModel and x a set of variables for edges
function edge_ub_callback()
   for (i,j) in E
     e = (i,j)
      if i != 0 && get_value(model.optimizer, x[e]) > 1.001
         println("Adding edge ub cut for e = ", e)
         add_dynamic_constr!(model.optimizer, [x[e]], [1.0], <=, 1.0, "edge_ub")
      end
   end 
end 
add_cut_callback!(model, edge_ub_callback, "edge_ub")
```
"""
function add_cut_callback!(user_model::VrpModel, callback::Any, constr_name::String)
   user_model.callbacks[constr_name] = callback
end

function add_callbacks_to_optimizer(optimizer::VrpOptimizer)

   cb_idx = 1
   for (constr_name, user_cb_fct) in optimizer.user_model.callbacks
      #auxiliary function
      aux_fct = function (cb)
         cb_info = optimizer.callbacks[constr_name]
         idx = getcstrindex(cb)[1] - cb_info.nb_added_constrs + length(cb_info.to_add_constrs)
         user_constr_info = cb_info.to_add_constrs[idx]
         for idx = 1:length(user_constr_info.vars)
            user_var, coeff = user_constr_info.vars[idx], user_constr_info.coeffs[idx]
            for var in optimizer.user_var_to_vars[user_var]
               addtermtoconstraint!(cb, var, coeff)
            end
         end
         addrhstoconstraint!(cb, user_constr_info.sense, user_constr_info.rhs)
      end
      if cb_idx == 1
         @dynconstraint(optimizer.formulation, bdone[], aux_fct)
      elseif cb_idx == 2
         @dynconstraint(optimizer.formulation, bdtwo[], aux_fct)
      elseif cb_idx == 3
         @dynconstraint(optimizer.formulation, bdthree[], aux_fct)
      elseif cb_idx == 4
         @dynconstraint(optimizer.formulation, bdfour[], aux_fct)
      elseif cb_idx == 5
         @dynconstraint(optimizer.formulation, bdfive[], aux_fct)
      elseif cb_idx == 6
         @dynconstraint(optimizer.formulation, bdsix[], aux_fct)
      elseif cb_idx == 7
         @dynconstraint(optimizer.formulation, bdseven[], aux_fct)
      elseif cb_idx == 8
         @dynconstraint(optimizer.formulation, bdeight[], aux_fct)
      elseif cb_idx == 9
         @dynconstraint(optimizer.formulation, bdnine[], aux_fct)
      elseif cb_idx == 10
         @dynconstraint(optimizer.formulation, bdten[], aux_fct)
      else
         error("VrpSolver supports at most 10 cut callbacks")
      end
      #separation function
      sep_fct = function (cb)
         #user callback adds cuts to cb_info
         user_cb_fct()
         cb_info = optimizer.callbacks[constr_name]
         for idx in 1:length(cb_info.to_add_constrs)
            if cb_info.cb_idx == 1
               @augmentconstraint(cb, bdone[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 2
               @augmentconstraint(cb, bdtwo[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 3
               @augmentconstraint(cb, bdthree[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 4
               @augmentconstraint(cb, bdfour[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 5
               @augmentconstraint(cb, bdfive[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 6
               @augmentconstraint(cb, bdsix[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 7
               @augmentconstraint(cb, bdseven[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 8
               @augmentconstraint(cb, bdeight[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 9
               @augmentconstraint(cb, bdnine[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            elseif cb_info.cb_idx == 10
               @augmentconstraint(cb, bdten[cb_info.nb_added_constrs - length(cb_info.to_add_constrs) + idx])
            else
	       error("VrpSolver supports at most 10 cut callbacks")
            end
         end
         empty!(cb_info.to_add_constrs)
      end

      addlazycallback(optimizer.formulation, sep_fct)
      optimizer.callbacks[constr_name] = CallbackInfo(constr_name, sep_fct, aux_fct, DynamicConstrInfo[], 0, cb_idx)
      cb_idx += 1
   end
end

"""
    add_capacity_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)

Define a *rounded capacity cut* (RCC) separator over a collection of packing sets defined on vertices.
RCC separators cannot be used if the packings sets are defined on arcs.

# Arguments
- `model::VrpModel`: model to be added the RCC separator.
- `demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}`: array of pairs of packing set and demand.
- `capacity::Float64`: capacity considered in the RCC separator.

# Examples
```julia
# Let PS be an array of `n` packing sets in vertices, where PS[i] is the i-th packing set
# Let d[i] be the demand associated with the i-th packing set
# Let `model` be a VrpModel and Q the capacity to be used in the RCC
add_capacity_cut_separator!(model, [(PS[i], d[i]) for i in 1:n], Q) # add a RCC separator
```
"""
function add_capacity_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)
   for (ps_set,d) in demands
      !(ps_set in model.packing_sets) && error("Collection that is not a packing set was used in a capacity cut separator." *
                                                " Only the packing set collections can be used for add_capacity_cut_separator")
   end

   # create and map variables to all uncovered arcs connecting packing set pairs
   if !isempty(model.arcs_by_packing_set_pairs)
      id_demands = [0 for i in 1:length(model.packing_sets)]
      for (ps_set, d) in demands
         ps_id = findall(x->x==ps_set, model.packing_sets)
         id_demands[ps_id[1]] = Int(d)
      end

      num_missing_arcs = 0
      uncovered = Tuple{Int, Int}[]
      dims_psp = size(model.arcs_by_packing_set_pairs)
      arcs_by_psp = model.arcs_by_packing_set_pairs
      for head in 1:dims_psp[1], tail in (head + 1):dims_psp[2]
         if (id_demands[head] > 0) && (id_demands[tail] > 0) && !isempty(arcs_by_psp[head,tail])
            push!(uncovered, (head, tail))
            num_missing_arcs += length(arcs_by_psp[head,tail])
         end
      end
      if length(uncovered) > 0
         println("VrpSolver: adding $(length(uncovered)) internal variables mapping to ",
            "$num_missing_arcs arcs for use by capacity cuts"
         )
      end
      if num_missing_arcs > 0
         @variable(model.formulation,
            RCCsepX[ps_pair in uncovered], Int
         )
         for (head, tail) in uncovered
            for (graph, arc) in arcs_by_psp[head, tail]
               add_arc_var_mapping!(graph, arc.id, RCCsepX[(head, tail)])
            end
            arcs_by_psp[head, tail] = Array{Tuple{VrpGraph, VrpArc}}(undef, 0)
         end
      end
   end

   push!(model.cap_cuts_info, CapacityCutInfo(demands, capacity))
end

function add_capacity_cut_separators_to_optimizer(optimizer::VrpOptimizer)
   user_model = optimizer.user_model
   for cap_cut_info in user_model.cap_cuts_info
      demands = [0 for i in 1:length(user_model.packing_sets)]
      for (ps_set, d) in cap_cut_info.demands
         ps_id = findall(x->x==ps_set, user_model.packing_sets)
         demands[ps_id[1]] = Int(d)
      end
      add_rcsp_capacity_cuts!(optimizer.formulation, Int(cap_cut_info.capacity), demands, is_facultative = false, root_priority_level = 3.0)
   end
end

"""
    add_strongkpath_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)

Define a *strong k-path* (SKP) separator over a collection of packing sets defined on vertices.
SKP separators cannot be used if the packings sets are defined on arcs.

# Arguments
- `model::VrpModel`: model to be added the RCC separator.
- `demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}`: array of pairs of packing set and demand.
- `capacity::Float64`: capacity considered in the RCC separator.

# Examples
```julia
# Let PS be an array of `n` packing sets in vertices, where PS[i] is the i-th packing set
# Let d[i] be the demand associated with the i-th packing set
# Let `model` be a VrpModel and Q the capacity to be used in the RCC
add_capacity_cut_separator!(model, [(PS[i], d[i]) for i in 1:n], Q) # add a RCC separator
```
"""
function add_strongkpath_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)
   for (ps_set,d) in demands
      !(ps_set in model.packing_sets) && error("Collection that is not a packing set was used in a strong k-path cut separator." *
                                                " Only the packing set collections can be used for add_strongkpath_cut_separator")
   end

   # create and map variables to all uncovered arcs connecting packing set pairs
   if !isempty(model.arcs_by_packing_set_pairs)
      id_demands = [0 for i in 1:length(model.packing_sets)]
      for (ps_set, d) in demands
         ps_id = findall(x->x==ps_set, model.packing_sets)
         id_demands[ps_id[1]] = Int(d)
      end

      num_missing_arcs = 0
      uncovered = Tuple{Int, Int}[]
      dims_psp = size(model.arcs_by_packing_set_pairs)
      arcs_by_psp = model.arcs_by_packing_set_pairs
      for head in 1:dims_psp[1], tail in (head + 1):dims_psp[2]
         if (id_demands[head] > 0) && (id_demands[tail] > 0) && !isempty(arcs_by_psp[head,tail])
            push!(uncovered, (head, tail))
            num_missing_arcs += length(arcs_by_psp[head,tail])
         end
      end
      if length(uncovered) > 0
         println("VrpSolver: adding $(length(uncovered)) internal variables mapping to ",
            "$num_missing_arcs arcs for use by strong k-path cuts"
         )
      end
      if num_missing_arcs > 0
         @variable(model.formulation,
            RCCsepX[ps_pair in uncovered], Int
         )
         for (head, tail) in uncovered
            for (graph, arc) in arcs_by_psp[head, tail]
               add_arc_var_mapping!(graph, arc.id, RCCsepX[(head, tail)])
            end
            arcs_by_psp[head, tail] = Array{Tuple{VrpGraph, VrpArc}}(undef, 0)
         end
      end
   end

   push!(model.strongkpath_cuts_info, CapacityCutInfo(demands, capacity))
end

function add_strongkpath_cut_separators_to_optimizer(optimizer::VrpOptimizer)
   user_model = optimizer.user_model
   for cap_cut_info in user_model.strongkpath_cuts_info
      demands = [0 for i in 1:length(user_model.packing_sets)]
      for (ps_set, d) in cap_cut_info.demands
         ps_id = findall(x->x==ps_set, user_model.packing_sets)
         demands[ps_id[1]] = Int(d)
      end
      add_rcsp_strongkpath_cuts!(optimizer.formulation, Int(cap_cut_info.capacity), demands, is_facultative = false, root_priority_level = 1.0)
   end
end

function get_mapped_containers_names(user_model::VrpModel)
   user_form = user_model.formulation
   user_var_to_graphs = extract_user_var_to_graphs(user_model)
   user_vars = [JuMP.Variable(user_form, i) for i in 1:user_form.numCols]

   names = String[]
   for user_var in user_vars
      (var_container_name, var_id) = split_var_name(user_var)
      if haskey(user_var_to_graphs, user_var)
         if !(var_container_name in names)
            push!(names, var_container_name) 
         end
      end
   end
   return names
end

function get_ignored_containers_names(user_model::VrpModel, mapped_names)
   user_form = user_model.formulation
   user_var_to_graphs = extract_user_var_to_graphs(user_model)
   user_vars = [JuMP.Variable(user_form, i) for i in 1:user_form.numCols]

   names = String[]
   for user_var in user_vars
      (var_container_name, var_id) = split_var_name(user_var)
      if !haskey(user_var_to_graphs, user_var) #unmapped var
         if (var_container_name in mapped_names) && !(var_container_name in names) 
            push!(names, var_container_name) 
         end
      end
   end
   return names
end

function build_optimizer_formulation(user_model::VrpModel)

   user_form = user_model.formulation
   user_var_to_graphs = extract_user_var_to_graphs(user_model)
   user_vars = [JuMP.Variable(user_form, i) for i in 1:user_form.numCols]
   mapped_names = get_mapped_containers_names(user_model)
   ignored_names = get_ignored_containers_names(user_model, mapped_names)

   formulation = VrpBlockModel()

   #creating variables
   var_container_name_to_ids = Dict{String,Array{Any}}()
   for user_var in user_vars
      (var_container_name, var_id) = split_var_name(user_var)
      if !haskey(user_var_to_graphs, user_var) #master variable
         if !(var_container_name in mapped_names) 
            ids = [0]
         else
            ids = []    
         end
      else
         ids = user_var_to_graphs[user_var] #subproblem variable
      end
      for id in ids
         if !haskey(var_container_name_to_ids, var_container_name)
            var_container_name_to_ids[var_container_name] = [(id, var_id)]
         else
            push!(var_container_name_to_ids[var_container_name], (id, var_id))
         end
      end
   end

   #workaround to use eval(parse()) for generating variables
   for (var_container_name, vars_ids) in var_container_name_to_ids
      vars_container = JuMP.JuMPArray(JuMP.Array{JuMP.variabletype(formulation)}(JuMP.undef, JuMP.length(vars_ids)), (vars_ids,))
      for a in vars_ids
         JuMP.coloncheck(a)
         vars_container[a] = JuMP.constructvariable!(formulation, msg -> error("Error when VrpSolver constructed variable :", msg), -Inf, Inf, :Default, JuMP.EMPTYSTRING, NaN)
      end
      JuMP.pushmeta!(vars_container, :model, formulation)
      JuMP.push!(formulation.dictList, vars_container)
      JuMP.registervar(formulation, Symbol(var_container_name), vars_container)
      JuMP.storecontainerdata(formulation, vars_container, Symbol(var_container_name), (vars_ids,), JuMP.IndexPair[JuMP.IndexPair(:a, :vars_ids)], Expr(:copyast, :((QuoteNode(:(()))))))
   end

   user_var_to_vars = extract_user_var_to_vars(user_form, formulation, user_var_to_graphs, mapped_names)

   #setting category, lb and ub of the variables
   for user_var in user_vars
      for var in user_var_to_vars[user_var]
         setcategory(var, getcategory(user_var))
         setlowerbound(var, getlowerbound(user_var))
         setupperbound(var, getupperbound(user_var))
      end
   end

   #creating objective function
   if getobjectivesense(user_form) == :Max
      error("VrpSolver does not handle maximization problems")
   end
   integer_objective = true
   coeffs, vars = [], []
   for (coeff,user_var) in linearterms(user_form.obj.aff)
      if getcategory(user_var) == :Cont || modf(coeff)[1] != 0.0
         integer_objective = false
      end
      for var in user_var_to_vars[user_var]
         push!(coeffs, coeff)
         push!(vars, var)
      end
   end

   JuMP.setobjective(formulation, :Min, AffExpr(vars, coeffs, user_form.obj.aff.constant))

   #creating constraints
   for constr in user_form.linconstr
      coeffs = []
      vars = []
      for (coeff,user_var) in linearterms(constr.terms)
         for var in user_var_to_vars[user_var]
            push!(coeffs, coeff)
            push!(vars, var)
         end
      end
      JuMP.addconstraint(formulation, LinearConstraint(AffExpr(vars, coeffs, 0.0), constr.lb, constr.ub))
   end

   return (formulation, user_var_to_vars, var_container_name_to_ids, integer_objective, mapped_names, ignored_names)
end

function set_branching_priorities_in_optimizer(user_model::VrpModel, formulation, var_container_name_to_ids, user_var_to_vars)

   var_container_name_to_probs = Dict{String,Array{Int}}()
   for var_container_name in keys(var_container_name_to_ids)
      var_container_name_to_probs[var_container_name] = []
   end
   for (var_container_name, ids) in var_container_name_to_ids
      for var_id in ids
         prob = var_id[1]
         if !(prob in var_container_name_to_probs[var_container_name])
            push!(var_container_name_to_probs[var_container_name], prob)
         end
      end
   end

   #no subproblem branching
   for var_container_name in keys(var_container_name_to_ids)
      for prob in var_container_name_to_probs[var_container_name]
         branchingpriorityinsubproblem(formulation[Symbol(var_container_name)], (prob > 0 ? :DW_SP : :DW_MASTER , prob), -1)
      end
   end

   #branching in master
   for var_container_name in keys(var_container_name_to_ids)
      probs = var_container_name_to_probs[var_container_name]
      if !haskey(user_model.branching_priorities, var_container_name)
         #no branching, priority = -1
         for prob in probs
            branchingpriorityinmaster(formulation[Symbol(var_container_name)], (prob > 0 ? :DW_SP : :DW_MASTER, prob), -1)
         end
      else
         priority =  user_model.branching_priorities[var_container_name]
         if 0 in probs
            branchingpriorityinmaster(formulation[Symbol(var_container_name)], (:DW_MASTER, 0), priority)
         end
         #now we decide if subproblem variables are aggregated or not
	      aggregate = (length(probs) >= 3) || (!(0 in probs) && (length(probs) >= 2))
         if aggregate
	         addbranching(formulation, :aggrsubprob_variables, Symbol(var_container_name),
		      priority = priority, highest_priority = 0.5, number_of_ignored_indices = 1, preprocessing = false)
            #if we are aggregating, we disable branching on a single variable
	         for prob in probs
	            if prob != 0
                  branchingpriorityinmaster(formulation[Symbol(var_container_name)], (:DW_SP, prob), -1)
	            end
	         end
	      else
	         for prob in probs
	            if prob != 0
                  branchingpriorityinmaster(formulation[Symbol(var_container_name)], (:DW_SP, prob), priority)
	            end
	         end
	      end
      end
   end

   #now we add branching on expressions
   #expression families
   for (exp_family, name, priority) in user_model.branching_exp_families
      formulation.ext[:branching_expression][Symbol(name)] = Array{Tuple{Array, Array, Array, Float64}, 1}() #user index, coeffs, vars, priority
      for index in keys(exp_family)
         varids, coeffs = [], Float64[]
         for user_var_idx in 1:length(exp_family[index...].vars)
            user_var = exp_family[index...].vars[user_var_idx]
            coeff = exp_family[index...].coeffs[user_var_idx]
            for var in user_var_to_vars[user_var]
               push!(varids, Cint(var.col))
               push!(coeffs, Cdouble(coeff))
            end
         end
         push!(formulation.ext[:branching_expression][Symbol(name)], (index, coeffs, varids, float(priority)))
      end
   end
   #single expressions
   for (exp, name, priority) in user_model.branching_exps
      formulation.ext[:branching_expression][Symbol(name)] = Array{Any, 1}()
      varids, coeffs = [], Float64[]
      for user_var_idx in 1:length(exp.vars)
         user_var = exp.vars[user_var_idx]
         coeff = exp.coeffs[user_var_idx]
         for var in user_var_to_vars[user_var]
            push!(varids, Cint(var.col))
            push!(coeffs, Cdouble(coeff))
         end
      end
      push!(formulation.ext[:branching_expression][Symbol(name)], ((1,), coeffs, varids, float(priority)))
   end
end

"""
    VrpOptimizer(user_model::VrpModel, param_file::String, instance_name = ""; baptreedot="BaPTree.dot")

Build an optimizer for a VrpModel.

# Arguments
- `model::VrpModel`: model to be considered
- `param_file::String`: path for the VRPSolver parameters file

# Optional arguments
- `instance_name::String`: the instance name to be shown in the results line (line with execution statistics).
- `baptreedot::String`: path to the file to output the BaP Tree in dot format.
"""
function VrpOptimizer(user_model::VrpModel, param_file::String, instance_name = ""; baptreedot="BaPTree.dot")

   #creating optimizer formulation
   formulation, user_var_to_vars, var_container_name_to_ids, integer_objective, 
   mapped_names, ignored_names = build_optimizer_formulation(user_model)

   #creating bapcod networks
   networks = generate_pricing_networks(user_model, user_var_to_vars)

   #dw configuration
   dw(cstrname, cstrmid) = (:DW_MASTER, 0)
   add_Dantzig_Wolfe_decomposition(formulation, dw)
   dw_vars(varname, varmid) = varmid[1][1] == 0 ? (:DW_MASTER, 0) : (:DW_SP, varmid[1][1])
   add_Dantzig_Wolfe_decomposition_on_variables(formulation, dw_vars)

   #registering rcsp solvers
   rcsp_oracle(subproblem) = networks[subproblem[3][1]]
   for graph in user_model.graphs
      generate_rcsp_oracle_for_dw_sp!(formulation, (:DW_SP, graph.id), rcsp_oracle)
   end

   #multiplicity of subproblems
   sp_mult(spid, sptype) = user_model.graphs[spid[1]].multiplicity
   addspmultiplicity(formulation, sp_mult)

   #branching priorities
   set_branching_priorities_in_optimizer(user_model, formulation, var_container_name_to_ids, user_var_to_vars)

   if user_model.use_rank1_cuts
      add_rcsp_lim_mem_rank1_cuts!(formulation)
   end
   if haskey(user_model.branching_priorities, "_res_cons_branching_")
      priority = user_model.branching_priorities["_res_cons_branching_"]
      add_rcsp_elemset_resource_consumption_branching!(formulation, priority)
   end
   if haskey(user_model.branching_priorities, "_ryanfoster_branching_")
      priority = user_model.branching_priorities["_ryanfoster_branching_"]
      add_rcsp_packset_ryanfoster_branching!(formulation, priority)
   end

   for name in ignored_names
      println("VRPSolver warning: unmapped $name vars were ignored because there are mapped $name vars")
      flush(stdout)
   end
   optimizer = VrpOptimizer(user_model, param_file, instance_name, 
	                    formulation, user_var_to_vars,
			    Dict{String,CallbackInfo}(),
			    Dict{JuMP.Variable,Float64}[],
			    integer_objective, -1,
                            mapped_names, ignored_names, Dict(), baptreedot)
   user_model.optimizer = optimizer

   add_callbacks_to_optimizer(optimizer)
   add_capacity_cut_separators_to_optimizer(optimizer)
   add_strongkpath_cut_separators_to_optimizer(optimizer)

   return optimizer
end

"""
    optimize!(optimizer::VrpOptimizer)

Solve a VRPSolver problem.

It returns a pair with the status of the execution and a flag indicating whether a solution was found, respectively.

# Example
```julia
# let model be a VrpModel
optimizer = VrpOptimizer(model, "path_to_config/config.cfg")
(status, solution_found) = optimize!(optimizer)
```
"""
function optimize!(optimizer::VrpOptimizer)
   optimizer.formulation.solver = BaPCodSolver(
	     param_file = optimizer.param_file,
	     integer_objective = optimizer.integer_objective, baptreedot_file = optimizer.baptreedot_file,
	     user_params = "--MaxNbOfStagesInColGenProcedure 3 --colGenSubProbSolMode 3 --MipSolverMultiThread 1 --ApplyStrongBranchingEvaluation true" )
   status = solve(optimizer.formulation)
   has_solution = register_solutions(optimizer)

   println("statistics_cols: instance & :Optimal & cutoff & :bcRecRootDb & :bcTimeRootEval & :bcCountNodeProc & :bcRecBestDb & :bcRecBestInc & :bcTimeMain \\\\")
   print("statistics: $(optimizer.instance_name) & ")
   print("$(status == :Optimal ? 1 : 0) & ")
   print("$(optimizer.initUB) & ")
   @printf("%.2f & ", getstatistic(optimizer.formulation, :bcRecRootDb))
   @printf("%.2f & ", getstatistic(optimizer.formulation, :bcTimeRootEval)/100)
   print("$(getstatistic(optimizer.formulation, :bcCountNodeProc)) & ")
   @printf("%.2f & ", getstatistic(optimizer.formulation, :bcRecBestDb))
   if has_solution
      if optimizer.integer_objective
         print("$(Int(floor(getstatistic(optimizer.formulation, :bcRecBestInc) + 0.5))) & ")
      else
         @printf("%.2f & ", getstatistic(optimizer.formulation, :bcRecBestInc))
      end
      optimizer.stats[:bcRecBestInc] = getstatistic(optimizer.formulation, :bcRecBestInc)
   else
      print("-- & ")
   end
   @printf("%.2f \\\\\n", getstatistic(optimizer.formulation, :bcTimeMain)/100)
   flush(stdout)

   optimizer.stats[:bcRecRootDb] = getstatistic(optimizer.formulation, :bcRecRootDb)
   optimizer.stats[:bcTimeRootEval] = getstatistic(optimizer.formulation, :bcTimeRootEval)
   optimizer.stats[:bcCountNodeProc] = getstatistic(optimizer.formulation, :bcCountNodeProc)
   optimizer.stats[:bcRecBestDb] = getstatistic(optimizer.formulation, :bcRecBestDb)
   optimizer.stats[:bcTimeMain] = getstatistic(optimizer.formulation, :bcTimeMain)
   optimizer.stats[:bcCountMastSol] = getstatistic(optimizer.formulation, :bcCountMastSol)
   optimizer.stats[:bcCountCol] = getstatistic(optimizer.formulation, :bcCountCol)
   optimizer.stats[:bcCountCutInMaster] = getstatistic(optimizer.formulation, :bcCountCutInMaster)
   optimizer.stats[:bcTimeMastMPsol] = getstatistic(optimizer.formulation, :bcTimeMastMPsol)
   optimizer.stats[:bcTimeCgSpOracle] = getstatistic(optimizer.formulation, :bcTimeCgSpOracle)
   optimizer.stats[:bcTimeCutSeparation] = getstatistic(optimizer.formulation, :bcTimeCutSeparation)
   optimizer.stats[:bcTimeAddCutToMaster] = getstatistic(optimizer.formulation, :bcTimeAddCutToMaster)
   optimizer.stats[:bcTimeSetMast] = getstatistic(optimizer.formulation, :bcTimeSetMast)
   optimizer.stats[:bcTimeRedCostFixAndEnum] = getstatistic(optimizer.formulation, :bcTimeRedCostFixAndEnum)
   optimizer.stats[:bcTimeEnumMPsol] = getstatistic(optimizer.formulation, :bcTimeEnumMPsol)
   optimizer.stats[:bcTimeSBphase1] = getstatistic(optimizer.formulation, :bcTimeSBphase1)
   optimizer.stats[:bcTimeSBphase2] = getstatistic(optimizer.formulation, :bcTimeSBphase2)
   optimizer.stats[:bcTimePrimalHeur] = getstatistic(optimizer.formulation, :bcTimePrimalHeur)

   return status, has_solution
end

function register_solutions(optimizer::VrpOptimizer)
   cols_in_sol = Tuple{Float64,Dict{JuMP.Variable,Float64}}[]
   bapcodsol = optimizer.formulation.internalModel.inner.solution

   # Initialize Columns
   BaPCod.getstartsolution(optimizer.formulation.internalModel.inner, bapcodsol)
   status = BaPCod.getnextsolution(bapcodsol) # ignore master solution
   while status == 1
	   m = BaPCod.getsolutionmultiplicity(bapcodsol)
	   push!(cols_in_sol, (m, Dict{JuMP.Variable,Float64}()))
	   status = BaPCod.getnextsolution(bapcodsol)
   end

   BaPCod.getstartsolution(optimizer.formulation.internalModel.inner, bapcodsol)
   status = BaPCod.getnextsolution(bapcodsol) # ignore master solution
   col_id = 1
   while status == 1
	   for (user_var, vars) in optimizer.user_var_to_vars
		   user_var_val = 0.0
		   for var in vars
			   value = Ref{Cdouble}(0.0)
			   BaPCod.c_getValueOfVar(optimizer.formulation.internalModel.inner.ptr, bapcodsol, Cint(var.col-1), value)
			   user_var_val += value[]
		   end
		   if user_var_val > 0
			   cols_in_sol[col_id][2][user_var] = user_var_val
		   end
	   end
	   status = BaPCod.getnextsolution(bapcodsol)
	   col_id += 1
   end

   optimizer.cols_in_sol = cols_in_sol
   if length(cols_in_sol) == 0
      return false
   else
      return true
   end
end

"""
    set_cutoff!(optimizer::VrpOptimizer, ub::Float64)

Set an upper bound (primal bound) for the execution.

"""
function set_cutoff!(optimizer::VrpOptimizer, ub::Float64)
    objectivevalueupperbound(optimizer.formulation, ub)
    optimizer.initUB = ub
    absolute_ub = abs(ub)
    if (absolute_ub < 10000.0)
        objectivevaluemagnitude(optimizer.formulation, 10000.0)
    else
        objectivevaluemagnitude(optimizer.formulation, absolute_ub)
    end
end
set_cutoff!(optimizer::VrpOptimizer, ub::Int) = set_cutoff!(optimizer, Float64(ub))

"""
    get_objective_value(optimizer::VrpOptimizer)

Get the objective function value after optimization.

"""
function get_objective_value(optimizer::VrpOptimizer)
   value = getobjectivevalue(optimizer.formulation)
   if optimizer.integer_objective
      return round(value)
   else
      return value
   end
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.Variable)

Get the value for a decision variable after optimization.

"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.Variable)
   val = 0.0
   for var in optimizer.user_var_to_vars[user_var]
      val += JuMP.getvalue(var)
   end
   return val
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.Variable})

Get the values for an array of decision variables after optimization.
   

"""
function get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.Variable})
   return [get_value(optimizer, user_var) for user_var in user_vars]
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.Variable, path_id::Int)

Get the value for a decision variable due to a specific path.

`path_id` shoul be a value between 1 and the number of positive paths.
"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.Variable, path_id::Int)
   if haskey(optimizer.cols_in_sol[path_id][2], user_var)
      return optimizer.cols_in_sol[path_id][2][user_var]
   else
      return 0.0
   end
end

"""
    get_value(optimizer::VrpOptimizer, path_id::Int)

Get the value for a path variable (λ variable created internally due to the mapping).

"""
function get_value(optimizer::VrpOptimizer, path_id::Int)
	if 1 <= path_id <= length(optimizer.cols_in_sol)
		return optimizer.cols_in_sol[path_id][1]
	end
	return 0.0
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.Variable}, path_id::Int)

Get the values for an array of decision variables that would be obtained from mapping only a single specific path variable (with value 1) to them. 
Necessary for identifying which paths are part of the solution in some models.

"""
function get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.Variable}, path_id::Int)
   return [get_value(optimizer, user_var, path_id) for user_var in user_vars]
end

"""
    get_number_of_positive_paths(optimizer::VrpOptimizer)

Get the number of paths (lambda variables) with positive value in the solution. Those paths will be numbered from 1 to get_number_of_positive_paths for the purpose 
of retrieving the value of the lambda variables and for identifying the path.

"""
function get_number_of_positive_paths(optimizer::VrpOptimizer)
   return length(optimizer.cols_in_sol)
end

"""
    add_dynamic_constr!(optimizer::VrpOptimizer, vars, coeffs, sense, rhs, constr_name::String)

Add user cut (dynamic constraint) to the formulation.

It should be called inside a cut callback function registered with [`add_cut_callback!`](@ref).

# Arguments
- `vars`: array of variables of the cut
- `coeffs`: array of coefficients, where `coeffs[i]` is the coefficient of `vars[i]`  
- `sense`: sense of the cut: <=, >= or ==. 
- `rhs`: right-hand side value
- `constr_name`: a name for the cut

# Example

```julia
# let model be a VrpModel and x a set of variables
function my_callback()
   # add cut 1.0*x[1] + 2.0*x[2] <= 1.0
   add_dynamic_constr!(model.optimizer, [x[1],x[2]], [1.0,2.0], <=, 1.0, "my_cut")
end 
add_cut_callback!(model, my_callback, "my_callback")
```
"""
function add_dynamic_constr!(optimizer::VrpOptimizer, vars, coeffs, sense, rhs, constr_name::String)
   constr_info = DynamicConstrInfo(vars, coeffs, sense, rhs)
   push!(optimizer.callbacks[constr_name].to_add_constrs, constr_info)
   optimizer.callbacks[constr_name].nb_added_constrs += 1
end

"""
    show(io::IO, graph::VrpGraph)

Show a VrpGraph.

"""
function show(io::IO, graph::VrpGraph)
   println(io, "===== Graph $(graph.id) =====")
   println(io, "L=$(graph.multiplicity[1]) U=$(graph.multiplicity[2])")
   println(io, "Resources")
   for r in graph.resources
      println(io, "$(r.id) $(r.is_main ? "main" : "secondary") $(r.is_disposable ? "disposable" : "non-disposable") $(r.is_binary ? "binary" : "")  $(r.is_automatic ? "automatic" : "")")
   end
   bin_res_legend = length(graph.vertices[1].res_bounds) > 0 ? "\nid (bin1, [lb_bin1,ub_bin1]) (bin2, [lb_bin2,ub_bin2]) ... (consumption interval for each binary resource)" : ""
   println(io, "Nodes$(bin_res_legend)")
   for v in graph.vertices
      if graph.cycle_problem && v.id == graph.sink_id
         continue
      end
      print(io, "$(v.user_id) ")
      for (i,rb) in enumerate(sort(collect(keys(v.res_bounds))))
         print(io, "($rb, [$(v.res_bounds[rb][1]),$(v.res_bounds[rb][2])]) ")
      end
      println(io)
   end
   println(io,"Arcs\n(i,j) arc_id consumption_res1 consumption_res2 ... (consumption interval for each non-binary resource) [list of mapped variables]")
   net_vertex_id_map = Dict(value => key for (key, value) in graph.user_vertex_id_map) # from user_id to net_id
   for a in graph.arcs
      i = net_vertex_id_map[a.tail]
      j = (graph.cycle_problem && a.head == graph.sink_id) ? net_vertex_id_map[graph.source_id] : net_vertex_id_map[a.head]
      print(io,"($i,$j) $(a.id) ")
      for rc in a.res_consumption
         print(io,"$rc ")
      end
      print(io,"(")
      for r in 1:length(a.res_consumption)
         if !haskey(a.res_bounds, r)
            print(io,"[,] ")
         else
            (lb, ub) = a.res_bounds[r]
            print(io,"[$lb, $ub] ")
         end
      end
      print(io,") [")
      for (id,mv) in enumerate(a.vars)
         coeff = mv[2] != 1 ? "$(mv[2])*" : ""
         if id < length(a.vars)
            print(io, string(coeff * "$(mv[1]),") )
         else
            print(io, string(coeff * "$(mv[1])") )
         end
      end
      print(io,"]\n")
   end
   println(io,"====================")
   flush(stdout)
end

function get_enum_paths(model::VrpModel, paramfile::String)
   optimizer = VrpOptimizer(model, paramfile, "")
   optimizer.formulation.ext[:complete_formulation] = true
   optimizer.formulation.solver = BaPCodSolver(
	      param_file = optimizer.param_file,
	      integer_objective = optimizer.integer_objective, baptreedot_file = optimizer.baptreedot_file,
	      user_params = "--MaxNbOfStagesInColGenProcedure 3 --colGenSubProbSolMode 3 --MipSolverMultiThread 1 --ApplyStrongBranchingEvaluation true" )
   solve(optimizer.formulation)
   paths = []
   bcsol = BaPCod.new!()
   nb = BaPCod.getenumeratedcols(optimizer.formulation.internalModel.inner, bcsol)
   status = BaPCod.getnextsolution(bcsol)
   while status == 1
      graph_id = BaPCod.c_get_prob_first_id(bcsol)
      bapcodarcids = BaPCod.c_get_sol_arcs_ids(bcsol)
      graph = model.graphs[graph_id]
      path = Int[]
      for bid in bapcodarcids
         rcsp_id = graph.net.bcid_2pos[bid] 
         arc_id = graph.arc_rcsp_id_to_id[rcsp_id] 
	      push!(path, arc_id)
	   end
	   push!(paths, (graph, path))
      status = BaPCod.getnextsolution(bcsol)
   end
   return paths
end

function print_enum_paths(model::VrpModel, paramfile::String)
    paths = get_enum_paths(model, paramfile)
    for (id, (graph, arcs)) in enumerate(paths)
        print("path $(id) graph $(graph.id): ")
         for arcid in arcs
	         print(arcid, " ") 
         end
	      println()
    end
end

"""
    print_enum_paths(paths)
   
Print the enumerated paths for all graphs.

Warning 1: the enumeration procedure only produces ``\\mathcal{E}``-elementarity-paths. Moreover, for all paths that visit exactly the same elementarity sets, only a single least cost path is produced.

# Example
```julia
# Let `model` be a VrpModel and `path_to_params_file` the path for the parameters file
enum_paths, complete_form = get_complete_formulation(model, path_to_params_file)
complete_form.solver = CplexSolver() # set MIP solver (it can be another one than CPLEX)
print_enum_paths(enum_paths)
```
"""
function print_enum_paths(paths)
   println("\n\nEnumerated paths (v1->(arc_id)->v2->...):")
   net_vertex_id_map = Dict{Int,Dict{Int,Int}}()
   for (id, (graph, arcs)) in enumerate(paths)
      if !haskey(net_vertex_id_map, graph.id)
         net_vertex_id_map[graph.id] = Dict(value => key for (key, value) in graph.user_vertex_id_map) # from user_id to net_id
      end
      print("path $(id) graph $(graph.id): ")
      arc = graph.arc_id_to_arc[arcs[1]]
      i = net_vertex_id_map[graph.id][arc.tail] 
      j = (graph.cycle_problem && arc.head == graph.sink_id) ? net_vertex_id_map[graph.id][graph.source_id] : net_vertex_id_map[graph.id][arc.head]
      print("$i-($(arcs[1]))->$j")
      prev = j 
      for arcid in arcs[2:end]
         arc = graph.arc_id_to_arc[arcid]
         j = (graph.cycle_problem && arc.head == graph.sink_id) ? net_vertex_id_map[graph.id][graph.source_id] : net_vertex_id_map[graph.id][arc.head]
         print("-($arcid)->$j") 
         prev = j
      end
      println()
   end
   println()
end

"""
    get_complete_formulation(model::VrpModel, paramfile::String)

Get the complete formulation, which includes mapping constraints with λ variables (for paths).

The enumerated paths and the complete formulation are returned. It has all enumerated paths for all graphs and a JuMP.Model. 
This function can be seen as a tool for debugging the VRPSolver model, for example, when applied to small instances to check the correctness
of the final model produced. The user may solve it with another solver option in JuMP 0.18 (e.g. for Gurobi, Coin CBC or CPLEX).

Warning 1: the enumeration procedure only produces ``\\mathcal{E}``-elementarity-paths. Moreover, for all paths that visit exactly the same elementarity sets, only a single least cost path is produced.

Warning 2: if some user cuts are essential for the correctness of the model (i.e., they are not only used to improve the linear relaxation), the corresponding cut callbacks should be called for solving the complete formulation.
    
# Example
```julia
# Let `model` be a VrpModel and `path_to_params_file` the path for the parameters file
enum_paths, complete_form = get_complete_formulation(model, path_to_params_file)
complete_form.solver = CplexSolver() # set MIP solver (it can be another one than CPLEX)
print_enum_paths(enum_paths)
println(complete_form)
solve(complete_form)
println("Objective value: ", getobjectivevalue(complete_form))
```
"""
function get_complete_formulation(model::VrpModel, paramfile::String)
   paths = get_enum_paths(model, paramfile)
   formulation = JuMP.Model()

   user_form = model.formulation
   user_var_to_graphs = extract_user_var_to_graphs(model)
   user_var_to_var = Dict{JuMP.Variable,JuMP.Variable}() 
   user_vars = [JuMP.Variable(user_form, i) for i in 1:user_form.numCols]
   mapped_names = get_mapped_containers_names(model)
   ignored_names = get_ignored_containers_names(model, mapped_names)

   function ignored_var(user_var)
      (var_container_name, var_id) = split_var_name(user_var)
      if (var_container_name in ignored_names) && !haskey(user_var_to_graphs, user_var)
          return true
      end
      return false
   end

   #creating variables
   vars_info = Dict{String,Array{Any}}()
   for user_var in user_vars
      if ignored_var(user_var)
         continue
      end
      (var_container_name, var_id) = split_var_name(user_var)
      if length(var_id) == 1
         var_id = var_id[1]
      end
      if !haskey(vars_info, var_container_name)
         vars_info[var_container_name] = []
      end
      push!(vars_info[var_container_name], (var_id, user_var, getcategory(user_var),
         getlowerbound(user_var), getupperbound(user_var))
      )
   end
   vars_info["λ"] = [(i, nothing, :Int, 0, Inf) for i in 1:length(paths)] #lambda variables

   lambda_container = nothing
   for (var_container_name, info) in vars_info
      vars_ids = [i[1] for i in info]
      vars_container = JuMP.JuMPArray(JuMP.Array{JuMP.variabletype(formulation)}(JuMP.undef, JuMP.length(vars_ids)), (vars_ids,))
      if var_container_name == "λ"
         lambda_container = vars_container
      end
      for (id, user_var, cat, lb, ub) in info
         JuMP.coloncheck(id)
         vars_container[id] = JuMP.constructvariable!(formulation, msg -> error("Error when VrpSolver constructed variable :", msg), -Inf, Inf, :Default, JuMP.EMPTYSTRING, NaN)
         setcategory(vars_container[id], cat)
         setlowerbound(vars_container[id], lb)
         setupperbound(vars_container[id], ub)
         if user_var != nothing
            user_var_to_var[user_var] = vars_container[id]
         end
      end
      JuMP.pushmeta!(vars_container, :model, formulation)
      JuMP.push!(formulation.dictList, vars_container)
      JuMP.registervar(formulation, Symbol(var_container_name), vars_container)
      JuMP.storecontainerdata(formulation, vars_container, Symbol(var_container_name), (vars_ids,), JuMP.IndexPair[JuMP.IndexPair(:a, :vars_ids)], Expr(:copyast, :((QuoteNode(:(()))))))
   end

   #computing uservar_map for paths
   paths_uservar_map = []
   for (graph, path) in paths
      push!(paths_uservar_map, get_path_uservar_map(path, graph))
   end 
   #creating objective function
   if getobjectivesense(user_form) == :Max
      error("VrpSolver does not handle maximization problems")
   end
   coeffs, vars = [], []
   for (coeff,user_var) in linearterms(user_form.obj.aff)
      if ignored_var(user_var)
         continue
      end
      push!(coeffs, coeff)
      push!(vars, user_var_to_var[user_var])
   end
   JuMP.setobjective(formulation, :Min, AffExpr(vars, coeffs, user_form.obj.aff.constant))

   #creating constraints
   for constr in user_form.linconstr
      coeffs = []
      vars = []
      for (coeff,user_var) in linearterms(constr.terms)
         if ignored_var(user_var)
            continue
         end
         push!(coeffs, coeff)
         push!(vars, user_var_to_var[user_var])
      end
      JuMP.addconstraint(formulation, LinearConstraint(AffExpr(vars, coeffs, 0.0), constr.lb, constr.ub))
   end

   #mapping constraints
   for user_var in user_vars
      if haskey(user_var_to_graphs, user_var)
         coeffs = [1]
         vars = [user_var_to_var[user_var]]
         #lambda vars
         for (id, (graph, path)) in enumerate(paths)
            if haskey(paths_uservar_map[id], user_var)
               push!(coeffs, -paths_uservar_map[id][user_var])
               push!(vars, lambda_container[id])
            end
         end
         JuMP.addconstraint(formulation, LinearConstraint(AffExpr(vars, coeffs, 0.0), 0.0, 0.0))
      end
   end
   
   #convexity constraint
   for graph in model.graphs
      coeffs = []
      vars = []
      #lambda vars
      for (id, (g, path)) in enumerate(paths)
         if g == graph
            push!(coeffs, 1)
            push!(vars, lambda_container[id])
         end
      end
      JuMP.addconstraint(formulation, LinearConstraint(AffExpr(vars, coeffs, 0.0), -Inf, graph.multiplicity[2]))
      JuMP.addconstraint(formulation, LinearConstraint(AffExpr(vars, coeffs, 0.0), graph.multiplicity[1], Inf))
   end

   return paths, formulation
end

function get_path_uservar_map(path::Vector{Int}, graph::VrpGraph)
   uservar_map = Dict{JuMP.Variable,Float64}()
   for arc_id in path
      arc = graph.arcs[arc_id]
      for (uservar, coef) in arc.vars
         if !haskey(uservar_map, uservar)
            uservar_map[uservar] = coef
	 else
            uservar_map[uservar] += coef
	 end
      end
   end
   return uservar_map
end

function get_path_coef_in_constr(path::Vector{Int}, graph::VrpGraph, uservar_map::Dict{JuMP.Variable,Float64}, constr)
   path_coef = 0.0 
   for (coef,user_var) in linearterms(constr.terms)
      if haskey(uservar_map, user_var)
         path_coef += coef*uservar_map[user_var]
      end
   end
   return path_coef
end

function get_path_coef_in_obj(path::Vector{Int}, graph::VrpGraph, uservar_map::Dict{JuMP.Variable,Float64}, obj_fct)
   path_coef = 0.0 
   for (coef,user_var) in linearterms(obj_fct.aff)
      if haskey(uservar_map, user_var)
         path_coef += coef*uservar_map[user_var]
      end
   end
   return path_coef
end
 
function VrpBlockModel(;solver = JuMP.UnsetSolver())
    m = JuMP.Model(solver = solver)
    JuMP.setsolvehook(m, vrp_hook)

    m.ext[:complete_formulation] = false
    m.ext[:already_called] = false
  
    # Block decomposition data
    m.ext[:block_decomposition] = BlockDecomposition.BlockDecompositionData()
    # Variables & constraints report
    m.ext[:varcstr_report] = BlockDecomposition.VarCstrReport()
    # Priorities & multiplicities
    m.ext[:sp_mult_fct] = nothing
    m.ext[:sp_prio_fct] = nothing
  
    # Storage (to list all subproblems)
    m.ext[:sp_list_dw] = Dict{Tuple, Integer}()
    m.ext[:sp_list_b] = Dict{Tuple, Integer}()
    m.ext[:sp_tab] = nothing
  
    # Data sent to the solver
    m.ext[:cstrs_decomposition_list] = nothing
    m.ext[:vars_decomposition_list] = nothing
    m.ext[:sp_mult_tab] = nothing
    m.ext[:sp_prio_tab] = nothing
    m.ext[:objective_data] = BlockDecomposition.ObjectiveData(NaN, -Inf, Inf)
  
    m.ext[:var_branch_prio_dict] = Dict{Tuple{Symbol, Tuple, Symbol}, Cdouble}() # (varname, sp, where) => (priority)
    m.ext[:branching_rules] = Dict{Symbol, Any}()
    m.ext[:branching_expression] = Dict{Symbol, Array{Tuple{Tuple, Array, Array, Float64}}}()
  
    # Callbacks
    m.ext[:oracles] = Array{Tuple{Tuple, Symbol, Function}}(undef, 0)
  
    m.ext[:generic_vars] = Dict{Symbol, Tuple{JuMP.Variable, Function}}()
    m.ext[:generic_cstrs] = Dict{Int, Tuple{JuMP.JuMP.ConstraintRef, String, Function}}()
    m.ext[:cstrs_preproc] = Dict{Tuple{Symbol, Tuple}, Bool}() # (cstrname, sp) => bool
  
    m.ext[:facultative_cuts_cb] = nothing
    m.ext[:core_cuts_cb] = nothing
  
    # Columns counter for generic variables & constraints
    m.ext[:colscounter] = 0
    m.ext[:rowscounter] = 0
    m
end

function vrp_hook(model;
     suppress_warnings=false,
     relaxation=false,
     kwargs...)

    if model.ext[:already_called]
        error("ERROR.")
    end

    model.ext[:already_called] = true
    # Step 1 : Create variables & constraints report
    BlockDecomposition.report_cstrs_and_vars!(model)
    BlockDecomposition.create_cstrs_vars_decomposition_list!(model)
    BlockDecomposition.create_sp_tab!(model)
    BlockDecomposition.create_sp_mult_tab!(model)
    BlockDecomposition.create_sp_prio_tab!(model)

    # Step 2 : Send decomposition (& others) data to the solver
    # Cstrs decomposition : mandatory
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_constrs_decomposition!, :cstrs_decomposition_list, false)
    # Vars decomposition : mandatory
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_vars_decomposition!, :vars_decomposition_list, false)
    # Subproblems
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_sp_ids!, :sp_tab, true)
    # Subproblems multiplicities
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_sp_mult!, :sp_mult_tab, false)
    # Subproblems priorities
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_sp_prio!, :sp_prio_tab, false)
    # Variables branching priority
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_var_branching_prio!, :var_branch_prio_dict, false)
    # Oracles
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_oracles!, :oracles, false)
    # Branching rules
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_branching_rules!, :branching_rules, false)
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_branching_exp!, :branching_expression, false)
    # Core & facultative cuts callbacks
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_corecuts_cb!, :core_cuts_cb, false)
    BlockDecomposition.send_to_solver!(model, BlockDecomposition.set_facultativecuts_cb!, :facultative_cuts_cb, false)

    if applicable(BlockDecomposition.send_extras!, model) # works with BlockDecompositionExtras
        BlockDecomposition.send_extras!(model)
    end

    # Objective bounds and magnitude
    obj = model.ext[:objective_data]
    if applicable(BlockDecomposition.set_objective_bounds_and_magnitude!, model.solver, obj.magnitude, obj.lb, obj.ub)
        BlockDecomposition.set_objective_bounds_and_magnitude!(model.solver, obj.magnitude, obj.lb, obj.ub)
    end

    model.ext[:colscounter] = model.numCols
    model.ext[:rowscounter] = length(model.ext[:cstrs_decomposition_list])

    #if applicable(defineannotations, model, model.ext[:vars_decomposition_list])
    #  defineannotations(model, model.ext[:vars_decomposition_list])
    #end

    # Step 3 : Build + solve
    if model.ext[:complete_formulation]
      println("\e[42m complete formulation \e[00m")
      JuMP.build(model)
      return :Unsolved
    else
      println("\e[41m solve \e[00m")
      a = JuMP.solve(model, suppress_warnings=suppress_warnings,
            ignore_solve_hook=true,
            relaxation=relaxation)
      return a
    end
end

end

