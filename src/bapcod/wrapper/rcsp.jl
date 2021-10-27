function wbcr_new(c_model::Ptr{Cvoid}, sp_bctype::Integer, sp_bcid::Array, nb_nodes::Integer, nb_es::Integer, nb_ps::Integer, nb_cs::Integer)
  ptr = @bcr_ccall("new", Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Ptr{Int}, Cint, Cint, Cint, Cint),
                                      c_model, Cint(sp_bctype), sp_bcid, Cint(nb_nodes), Cint(nb_es), Cint(nb_ps), Cint(nb_cs))
end

function wbcr_delete(c_net::Ptr{Cvoid})
   @bcr_ccall("delete", Cvoid, (Ptr{Cvoid},), c_net)
end

function wbcr_new_resource(c_net::Ptr{Cvoid}, res_id::Integer)
  status = @bcr_ccall("newResource", Cint, (Ptr{Cvoid}, Cint),
                                            c_net, Cint(res_id))
  (status != 1) && error("Cannot create resource $res_id.")
end

function wbcr_set_as_main_resource(c_net::Ptr{Cvoid}, res_id::Integer, stepvalue::Cdouble)
  @bcr_ccall("setAsMainResource", Cvoid, (Ptr{Cvoid}, Cint, Cdouble),
                                         c_net, Cint(res_id), stepvalue)
end

function wbcr_set_as_nondisposable_resource(c_net::Ptr{Cvoid}, res_id::Integer)
  @bcr_ccall("setAsNonDisposableResource", Cvoid, (Ptr{Cvoid}, Cint),
                                         c_net, Cint(res_id))
end

function wbcr_set_special_as_nondisposable_resource(c_net::Ptr{Cvoid}, res_id::Integer)
  @bcr_ccall("setSpecialResourceAsNonDisposable", Cvoid, (Ptr{Cvoid}, Cint),
                                         c_net, Cint(res_id))
end

function wbcr_set_vertex_consumption_lb(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, lb::Cdouble)
  @bcr_ccall("setVertexConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                              c_net, Cint(n_id), Cint(res_id), lb)
end

function wbcr_set_vertex_consumption_ub(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, ub::Cdouble)
  @bcr_ccall("setVertexConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                              c_net, Cint(n_id), Cint(res_id), ub)
end

function wbcr_set_vertex_special_consumption_lb(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, lb::Cdouble)
  @bcr_ccall("setVertexSpecialConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                              c_net, Cint(n_id), Cint(res_id), lb)
end

function wbcr_set_vertex_special_consumption_ub(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, ub::Cdouble)
  @bcr_ccall("setVertexSpecialConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                              c_net, Cint(n_id), Cint(res_id), ub)
end

function wbcr_set_source(c_net::Ptr{Cvoid}, n_id::Integer)
  @bcr_ccall("setSource", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_set_sink(c_net::Ptr{Cvoid}, n_id::Integer)
  @bcr_ccall("setSink", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_new_arc(c_net::Ptr{Cvoid}, src::Integer, dst::Integer, cost::Cdouble)
  edge_id = @bcr_ccall("newArc", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                        c_net, Cint(src), Cint(dst), cost)
  return edge_id
end

function wbcr_attach_bcvar_to_arc(c_net::Ptr{Cvoid}, edge_id::Integer, c_model::Ptr{Cvoid}, var_col::Integer, var_coeff::Cdouble)
  @bcr_ccall("attachBcVarToArc", Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Cdouble),
                                        c_net, Cint(edge_id), c_model, Cint(var_col), var_coeff)
end

function wbcr_set_edge_consumption_lb(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Cdouble)
  @bcr_ccall("setArcConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                               c_net, Cint(edge_id), Cint(res_id), value)
end

function wbcr_set_edge_consumption_ub(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Cdouble)
  @bcr_ccall("setArcConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                               c_net, Cint(edge_id), Cint(res_id), value)
end

function wbcr_set_edge_consumption_value(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Cdouble)
  @bcr_ccall("setEdgeConsumptionValue", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                               c_net, Cint(edge_id), Cint(res_id), value)
end

function wbcr_set_edge_special_consumption_value(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Cdouble)
  @bcr_ccall("setEdgeSpecialConsumptionValue", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
                                               c_net, Cint(edge_id), Cint(res_id), value)
end

function wbcr_attach_elementarity_set_to_node(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
  @bcr_ccall("attachElementaritySetToNode", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                   c_net, Cint(n_id), Cint(es_id))
end

function wbcr_attach_elementarity_set_to_edge(c_net::Ptr{Cvoid}, edge_id::Integer, es_id::Integer)
  @bcr_ccall("attachElementaritySetToEdge", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                   c_net, Cint(edge_id), Cint(es_id))
end

function wbcr_add_vertex_to_mem_of_elementarity_set(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
  status = @bcr_ccall("addVertexToMemOfElementaritySet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                       c_net, Cint(n_id), Cint(es_id))
  (status != 1) && error("Cannot add vertex $n_id to memory of elementarity set $n_es.")
end

function wbcr_add_edge_to_mem_of_elementarity_set(c_net::Ptr{Cvoid}, edge_id::Integer, es_id::Integer)
  status = @bcr_ccall("addEdgeToMemOfElementaritySet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                      c_net, Cint(edge_id), Cint(es_id))
  (status != 1) && error("Cannot add edge $edge_id to memory of elementarity set $edge_id.")
end

function wbcr_add_vertex_to_packing_set(c_net::Ptr{Cvoid}, n_id::Integer, ps_id::Integer)
  status = @bcr_ccall("addVertexToPackingSet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                      c_net, Cint(n_id), Cint(ps_id))
end

function wbcr_add_edge_to_packing_set(c_net::Ptr{Cvoid}, edge_id::Integer, ps_id::Integer)
  status = @bcr_ccall("addEdgeToPackingSet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                      c_net, Cint(edge_id), Cint(ps_id))
end

function wbcr_add_vertex_to_covering_set(c_net::Ptr{Cvoid}, n_id::Integer, cs_id::Integer)
  status = @bcr_ccall("addVertexToCoveringSet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                      c_net, Cint(n_id), Cint(cs_id))
end

function wbcr_add_edge_to_covering_set(c_net::Ptr{Cvoid}, edge_id::Integer, cs_id::Integer)
  status = @bcr_ccall("addVertexToCoveringSet", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                      c_net, Cint(edge_id), Cint(cs_id))
end

function wbcr_add_packing_set_to_packing_set_cut_neighbourhood(c_net::Ptr{Cvoid}, ps2add::Integer, ps::Integer)
  status = @bcr_ccall("addPackingSetToPackingSetCutNeighbourhood", Cint, (Ptr{Cvoid}, Cint, Cint),
                                                                  c_net, Cint(ps2add), Cint(ps))
end

function wbcr_set_elementarity_sets_distance_matrix(c_net::Ptr{Cvoid}, dists::Array{Array{Cdouble,1},1}, nb_packsets::Integer)
  status = @bcr_ccall("setElemSetsDistanceMatrix", Cint, (Ptr{Cvoid}, Ptr{Ptr{Cdouble}}, Cint),
                                                                  c_net, dists, Cint(nb_packsets))
end

# function wbcr_add_elementarity_set_to_mem_of_elementarity_set_cut(c_net::Ptr{Cvoid}, es_id_2_add::Integer, es_id::Integer)
#   status = @bcr_ccall("addElementaritySetToMemOfElementaritySetCut", Cint, (Ptr{Cvoid}, Cint, Cint),
#                                                       c_net, Cint(es_id_2_add), Cint(es_id))
#   (status != 1) && error("Cannot add elem. set $es_id_2_add to memory of elem. set $es_id.")
# end

function wbcr_create_oracle(c_net::Ptr{Cvoid}, c_model::Ptr{Cvoid}, sp_type::Integer, sp_id::Array, save_standalone::Bool, standalone_filename::String)
  @bcr_ccall("createOracle", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cint}, UInt8, Ptr{UInt8}),
                                    c_net, c_model, sp_type, sp_id, save_standalone, standalone_filename)
end

function wbcr_add_generic_capacity_cut(c_model::Ptr{Cvoid}, max_cap::Integer,
  dem::Array, dem_length::Integer, is_fac::Integer, root_prio_lvl::Cdouble,
  non_root_prio_lvl::Cdouble)
  status = @bcr_ccall("addGenericCapacityCut", Cint, (Ptr{Cvoid}, Cint,
                                                      Ptr{Cint}, Cint, Cint, Cdouble, Cdouble),
                          c_model, Cint(max_cap), dem, Cint(dem_length), Cint(is_fac),
                          root_prio_lvl, non_root_prio_lvl)
  (status != 1) && error("Cannot add the generic capacity cut generator.")
end

function wbcr_add_generic_strongkpath_cut(c_model::Ptr{Cvoid}, max_cap::Integer,
  dem::Array, dem_length::Integer, is_fac::Integer, root_prio_lvl::Cdouble,
  non_root_prio_lvl::Cdouble)
  status = @bcr_ccall("addGenericStrongKPathCut", Cint, (Ptr{Cvoid}, Cint,
                                                      Ptr{Cint}, Cint, Cint, Cdouble, Cdouble),
                          c_model, Cint(max_cap), dem, Cint(dem_length), Cint(is_fac),
                          root_prio_lvl, non_root_prio_lvl)
  (status != 1) && error("Cannot add the generic strong k-path cut generator.")
end

function wbc_add_generic_clique_cut(c_model::Ptr{Cvoid})
  status = @bcr_ccall("addGenericCliqueCut", Cint, (Ptr{Cvoid},),
                                                        c_model)
  (status != 1) && error("Cannot add the generic clique cut.")
end

function wbc_add_generic_lim_mem_one_cut(c_model::Ptr{Cvoid})
  status = @bcr_ccall("addGenericLimMemOneCut", Cint, (Ptr{Cvoid},),
                                                        c_model)
  (status != 1) && error("Cannot add the generic lim-mem-one-rank cut.")
end

function wbcr_add_paths_per_network_branching(c_model::Ptr{Cvoid}, priority::Cdouble)
  status = @bcr_ccall("addPathsPerNetworkBranching", Cint, (Ptr{Cvoid}, Cdouble),
                                                            c_model, priority)
  (status != 1) && error("Cannot add the paths_per_network_branching.")
end

function wbcr_add_elemset_assign_branching(c_model::Ptr{Cvoid}, priority::Cdouble)
  status = @bcr_ccall("addElemSetAssignBranching", Cint, (Ptr{Cvoid}, Cdouble),
                                                          c_model, priority)
  (status != 1) && error("Cannot add the elemset_assign_branching.")
end

function wbcr_add_packset_ryanfoster_branching(c_model::Ptr{Cvoid}, priority::Cdouble)
  status = @bcr_ccall("addPackSetRyanAndFosterBranching", Cint, (Ptr{Cvoid}, Cdouble),
                                                          c_model, priority)
  (status != 1) && error("Cannot add the packset_ryanfoster_branching.")
end

function wbcr_add_permanent_ryanfoster_constraint(c_net::Ptr{Cvoid}, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool)
  status = @bcr_ccall("addPermanentRyanAndFosterConstraint", Cint, (Ptr{Cvoid}, Cint, Cint, UInt8),
                                                                  c_net, firstPackSetId, secondPackSetId, together)
end

function wbcr_add_elemset_resource_consumption_branching(c_model::Ptr{Cvoid}, priority::Cdouble)
  status = @bcr_ccall("addElemSetResourceConsumptionBranching", Cint, (Ptr{Cvoid}, Cdouble),
                                                          c_model, priority)
  (status != 1) && error("Cannot add the elemset_resource_consumption_branching.")
end

function wbcr_add_associated_var_to_resource(c_net::Ptr{Cvoid}, res_id::Integer, c_model::Ptr{Cvoid}, col_id::Integer)
  status = @bcr_ccall("addAssociatedVarToResource", Cint,
                      (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint),
                      c_net, Cint(res_id), c_model, Cint(col_id))
  (status != 1) && error("Cannot associate the variable with col $col_id to resource with id $res_id.")
end 

function wbcr_add_generic_knapsack_consumption_cut(c_model::Ptr{Cvoid}, rootpriolvl::Cdouble, nonrootpriolvl::Cdouble)
  status = @bcr_ccall("addGenericKnapsackConsumptionCut", Cint,
                      (Ptr{Cvoid}, Cdouble, Cdouble),
                      c_model, rootpriolvl, nonrootpriolvl)
  (status != 1) && error("Cannot add generic knapsack consumption cut.")
end

function wbcr_add_special_resource_consumption(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, consumption::Integer, lb::Integer, ub::Integer)
  status = @bcr_ccall("addSpecialResourceConsumption", Cint,
                      (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint),
                      c_net, Cint(n_id), Cint(res_id), Cint(consumption), Cint(lb), Cint(ub))
end
