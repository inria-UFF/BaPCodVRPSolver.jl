function bj_solve(model;
                suppress_warnings=false,
                relaxation=false,
                kwargs...)
  # Step 1 : Create variables & constraints report
  report_cstrs_and_vars!(model)
  create_cstrs_vars_decomposition_list!(model)
  create_sp_tab!(model)
  create_sp_mult_tab!(model)
  create_sp_prio_tab!(model)

  # Step 2 : Send decomposition (& others) data to the solver
  # Cstrs decomposition : mandatory
  send_to_solver!(model, set_constrs_decomposition!, :cstrs_decomposition_list, false)
  # Vars decomposition : mandatory
  send_to_solver!(model, set_vars_decomposition!, :vars_decomposition_list, false)
  # Subproblems
  send_to_solver!(model, set_sp_ids!, :sp_tab, true)
  # Subproblems multiplicities
  send_to_solver!(model, set_sp_mult!, :sp_mult_tab, false)
  # Subproblems priorities
  send_to_solver!(model, set_sp_prio!, :sp_prio_tab, false)
  # Variables branching priority
  send_to_solver!(model, set_var_branching_prio!, :var_branch_prio_dict, false)
  # Oracles
  send_to_solver!(model, set_oracles!, :oracles, false)
  # Branching rules
  send_to_solver!(model, set_branching_rules!, :branching_rules, false)
  send_to_solver!(model, set_branching_exp!, :branching_expression, false)
  # Core & facultative cuts callbacks
  send_to_solver!(model, set_corecuts_cb!, :core_cuts_cb, false)
  send_to_solver!(model, set_facultativecuts_cb!, :facultative_cuts_cb, false)

  if applicable(send_extras!, model) # works with BlockDecompositionExtras
    send_extras!(model)
  end

  # Objective bounds and magnitude
  obj = model.ext[:objective_data]
  if applicable(set_objective_bounds_and_magnitude!, model.solver, obj.magnitude, obj.lb, obj.ub)
    set_objective_bounds_and_magnitude!(model.solver, obj.magnitude, obj.lb, obj.ub)
  end

  model.ext[:colscounter] = model.numCols
  model.ext[:rowscounter] = length(model.ext[:cstrs_decomposition_list])

  #if applicable(defineannotations, model, model.ext[:vars_decomposition_list])
  #  defineannotations(model, model.ext[:vars_decomposition_list])
  #end

  # Step 3 : Build + solve
  a = JuMP.solve(model, suppress_warnings=suppress_warnings,
                    ignore_solve_hook=true,
                    relaxation=relaxation)
  a
end

function create_sp_tab!(model::JuMP.Model)
  for (sp_id, sp_type, f) in model.ext[:oracles]
    list_sp!(model, sp_type, sp_id)
  end
  model.ext[:sp_tab] = Vector{Tuple{Tuple, Symbol}}()
  for spid in collect(keys(model.ext[:sp_list_dw]))
    push!(model.ext[:sp_tab], (spid, :DW_SP))
  end
  for spid in collect(keys(model.ext[:sp_list_b]))
    push!(model.ext[:sp_tab], (spid, :B_SP))
  end
end

function send_to_solver!(model::JuMP.Model, f::Function, k::Symbol, mandatory::Bool)
  if haskey(model.ext, k)
    if applicable(f, model.solver, model.ext[k])
      f(model.solver, model.ext[k])
    else
      if mandatory
        @warn("Your solver does not support function $f.")
      end
    end
  else
    bjerror("Key $(k) does not exist in model.ext.")
  end
end

## Send data to BlockDecompositionExtras package
send_extras!() = nothing

function send_extras_to_solver!(model) end
export send_extras_to_solver!

# @require BlockDecompositionExtras begin
  function send_extras!(model::JuMP.Model)
    # if applicable(BlockDecompositionExtras.send_extras_to_solver!, model)
    #   BlockDecompositionExtras.send_extras_to_solver!(model)
    # end
    send_extras_to_solver!(model)
  end
# end

getdisaggregatedvalue(x::JuMP.JuMPContainer) = @warn("getdisaggregatedvalue of a JuMPContainer no longer available. Use getdisaggregatedvalue(x::JuMP.Variable).")
getdisaggregatedvalue(model::JuMP.Model) = @warn("getdisaggregatedvalue of a JuMP.Model no longer available. Use getdisaggregatedvalue(x::JuMP.Variable).")

function getdisaggregatedvalue(x::JuMP.Variable)
  if applicable(getdisaggregatedvalueofvariable, x.m.internalModel, x.col)
    return getdisaggregatedvalueofvariable(x.m.internalModel, x.col)
  else
    @warn("Your solver seems to not support disaggregated solutions.")
  end
end
