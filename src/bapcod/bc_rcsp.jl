function register_rcsp_oracles!(m)
  sp_bcid = Array{Cint}(undef, 8)
  for ((sp_type, sp_id), rcsp_fct) in m.solver.rcsp_oracles
    (sp_type != :DW_SP) && error("RCSP oracles only work on Dantzig-Wolfe subproblems.")

    network = rcsp_fct((m.inner.ptr, sp_type, sp_id))

    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    register_network!(m, sp_bctype, sp_bcid, network)
    wbcr_create_oracle(network.c_net.ptr, m.inner.ptr, sp_bctype, sp_bcid, network.save_standalone, network.standalone_filename)
  end
end

mutable struct C_network
  ptr::Ptr{Cvoid}
  function C_network(nb_nodes::Integer, nb_es::Integer, nb_ps::Integer, nb_cs::Integer, c_model::Ptr{Cvoid}, sp_bctype, sp_bcid)
    ptr = wbcr_new(c_model, sp_bctype, sp_bcid, nb_nodes, nb_es, nb_ps, nb_cs)
    c_net = new(ptr)
    finalizer(delete_C_network, c_net)
    return c_net
  end
end

function delete_C_network(c_net::C_network)
  if c_net.ptr != C_NULL
    wbcr_delete(c_net.ptr)
    c_net.ptr = C_NULL
  end
end

function register_network!(m, sp_bctype, sp_bcid, network)
  nb_resources = length(network.resources)
  nb_nodes = network.nb_nodes
  nb_es = network.nb_elementaritysets
  nb_ps = network.nb_packingsets
  nb_cs = network.nb_coveringsets
  # Step 1 : create the network object
  network.c_net = C_network(nb_nodes, nb_es, nb_ps, nb_cs, m.inner.ptr, sp_bctype, sp_bcid)
  c_net_ptr = network.c_net.ptr
  # Step 2 : register resources
  for (res_id, (main, disposable, step)) in network.resources
    wbcr_new_resource(c_net_ptr, res_id)
    main && wbcr_set_as_main_resource(c_net_ptr, res_id, step)
    if !disposable 
     wbcr_set_as_nondisposable_resource(c_net_ptr, res_id)
    end
  end
  #special resources
  for res_id in 1:length(network.specialresources)
     if !network.specialresources[res_id]
        wbcr_set_special_as_nondisposable_resource(c_net_ptr, res_id - 1)
     end
  end

  ## and node consumptions
  for n in 1:nb_nodes, (r, (lb, ub)) in network.verticesconsumptionsbounds[n]
    (lb != NaN) && wbcr_set_vertex_consumption_lb(c_net_ptr, n-1, r, lb)
    (ub != NaN) && wbcr_set_vertex_consumption_ub(c_net_ptr, n-1, r, ub)
  end
  ## and node special consumptions
  for n in 1:nb_nodes, (r, (lb, ub)) in network.verticesspecialconsumptionsbounds[n]
    (lb != NaN) && wbcr_set_vertex_special_consumption_lb(c_net_ptr, n-1, r-1, lb)
    (ub != NaN) && wbcr_set_vertex_special_consumption_ub(c_net_ptr, n-1, r-1, ub)
  end
  
  # Step 3 : register source & sink
  (network.source != 0) && wbcr_set_source(c_net_ptr, network.source - 1)
  (network.sink != 0) && wbcr_set_sink(c_net_ptr, network.sink - 1)
  # Step 4 : create edges and register attached variables
  edges_pos2bcid = zeros(Integer, length(network.edge_id2pos))
  for (pos, edge) in enumerate(network.edge_pos2id)
    s, d, o = edge
    # Step a : register the arc
    cost = Cdouble(network.edgescosts[pos])
    edge_id = wbcr_new_arc(c_net_ptr, s - 1, d - 1, cost)
    edges_pos2bcid[pos] = edge_id
    network.bcid_2pos[edge_id] = pos
    # Step b : register the variable
    variables = network.edgesvars[pos]
    if length(variables) != 0
      for var  in variables
        varcol = var[1] - 1
        varcoeff = Cdouble(var[2])
        wbcr_attach_bcvar_to_arc(c_net_ptr, edge_id, m.inner.ptr, varcol, varcoeff)
      end
    end
    # Step c : edge resources consumptions
    consumptions = network.edgesconsumptions[pos]
    for (r, value) in consumptions
      wbcr_set_edge_consumption_value(c_net_ptr, edge_id, r, Cdouble(value))
    end
    # Step d: edge special resource consumptions
    specialconsumptions = network.edgesspecialconsumptions[pos]
    for (r, value) in specialconsumptions
      wbcr_set_edge_special_consumption_value(c_net_ptr, edge_id, r-1, Cdouble(value))
    end
    # Step e: edge consumption bounds
    consumptionsbounds = network.arcsconsumptionsbounds[pos]
    for (r, (lb, ub)) in consumptionsbounds
      wbcr_set_edge_consumption_lb(c_net_ptr, edge_id, r, Cdouble(lb))
      wbcr_set_edge_consumption_ub(c_net_ptr, edge_id, r, Cdouble(ub))
    end
  end
  # Step 5 : create elementarity sets
  for es_id in 1:nb_es
    # attached to node
    for n_id in network.elemsets_ni_nodes[es_id]
      wbcr_attach_elementarity_set_to_node(c_net_ptr, n_id - 1, es_id - 1)
    end
    for edge_pos in network.elemsets_ni_edges[es_id]
      wbcr_attach_elementarity_set_to_edge(c_net_ptr, edges_pos2bcid[edge_pos], es_id - 1)
    end
  end
  # Step 6 : mem of elem sets
  for (es_pos, vertices) in enumerate(network.moes_vertex), n in vertices
    wbcr_add_vertex_to_mem_of_elementarity_set(c_net_ptr, n - 1, es_pos - 1)
  end
  for (es_pos, edges) in enumerate(network.moes_edges), e in edges
    wbcr_add_edge_to_mem_of_elementarity_set(c_net_ptr, edges_pos2bcid[e], es_pos - 1)
  end
  # Step 7 : packing sets
  for ps_id in 1:nb_ps
    for n_id in network.packsets_ni_nodes[ps_id]
      wbcr_add_vertex_to_packing_set(c_net_ptr, n_id - 1, ps_id - 1)
    end
    for edge_pos in network.packsets_ni_edges[ps_id]
      wbcr_add_edge_to_packing_set(c_net_ptr, edges_pos2bcid[edge_pos], ps_id - 1)
    end
  end
  # Step 8 : packing sets cuts neighbourhoods
  for (ps_id, psets2add) in enumerate(network.ps2pscutneighbourhood), ps2add in psets2add
    wbcr_add_packing_set_to_packing_set_cut_neighbourhood(c_net_ptr, ps2add - 1, ps_id - 1)
  end
  # Step 9 : covering sets
  for cs_id in 1:nb_cs
    for n_id in network.covsets_ni_nodes[cs_id]
      wbcr_add_vertex_to_covering_set(c_net_ptr, n_id - 1, cs_id - 1)
    end
    for edge_pos in network.covsets_ni_edges[cs_id]
      wbcr_add_edge_to_covering_set(c_net_ptr, edges_pos2bcid[edge_pos], cs_id - 1)
    end
  end
  #Step 10: distance matrix
  if network.es_dist_matrix != nothing
     wbcr_set_elementarity_sets_distance_matrix(c_net_ptr, network.es_dist_matrix, nb_es)
  end
  # Step 11 : Variables associated to resource
  for (res_id, col_id) in network.resource_associated_var
    wbcr_add_associated_var_to_resource(c_net_ptr, res_id, m.inner.ptr, col_id)
  end
  #Step 12: Ryan and Foster permanent constraints
  for (ps1,ps2,tog) in network.ryanfoster_constraints
    wbcr_add_permanent_ryanfoster_constraint(c_net_ptr, ps1, ps2, tog)
  end
  # for (es_pos, es2add) in enumerate(network.moes_elemset), espos2add in es2add
  #   wbcr_add_elementarity_set_to_mem_of_elementarity_set_cut(c_net_ptr, espos2add - 1, es_pos - 1)
  # end
end

function register_rcsp_genericcuts!(m)
  sp_bcid = Array{Cint}(undef, 8)
  for cut in m.solver.rcsp_generic_cuts
    c_model = m.inner.ptr
    if cut[1] == :capacity
      (max_cap, dem, isfac, rpl, nrpl, tpcri) = cut[2]
      wbcr_add_generic_capacity_cut(c_model, max_cap, dem, length(dem), isfac,
                                    rpl, nrpl, tpcri)
    elseif cut[1] == :strongkpath                                    
      (max_cap, dem, isfac, rpl, nrpl) = cut[2]
      wbcr_add_generic_strongkpath_cut(c_model, max_cap, dem, length(dem), isfac,
                                       rpl, nrpl)
    elseif cut[1] == :limmemrank1
      wbc_add_generic_lim_mem_one_cut(c_model)
    elseif cut[1] == :clique
      wbc_add_generic_clique_cut(c_model)
    else
      @warn("Unrecognized RCSP generic cut type.")
    end
  end
end

function register_rcsp_branching!(m)
  c_model = m.inner.ptr
  for branching in m.solver.rcsp_branching
    if branching[1] == :elemset_assign_branching
      priority = Cdouble(branching[2])
      wbcr_add_elemset_assign_branching(c_model, priority)
    elseif branching[1] == :elemset_rescons_branching
      priority = Cdouble(branching[2])
      wbcr_add_elemset_resource_consumption_branching(c_model, priority)
    elseif branching[1] == :packset_ryanfoster_branching
      priority = Cdouble(branching[2])
      wbcr_add_packset_ryanfoster_branching(c_model, priority)
    elseif branching[1] == :paths_per_network_branching
      priority = Cdouble(branching[2])
      wbcr_add_paths_per_network_branching(c_model, priority)
    elseif branching[1] == :knpconsumption_branching
      rootprio, nonrootprio = branching[2]
      wbcr_add_generic_knapsack_consumption_cut(c_model, Cdouble(rootprio), Cdouble(nonrootprio))
    else
      @warn("Unrecognized RCSP branching type.")
    end
  end
end
