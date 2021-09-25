# Entry method for variables registration with variables names and indexes
function c_register_vars(m::BcModel, A, l, u, c, vars_decomposition)
  var_bcid = Array{Cint}(undef, 8)
  sp_bcid = Array{Cint}(undef, 8)
  rows = rowvals(A)

  for (column_id, (name, v_id, sp_type, sp_id)) in enumerate(vars_decomposition)
     # BaPCod needs an index
    v_id = (v_id == nothing) ? (0,) : v_id
    sp_id = (sp_id == nothing) ? (0,) : sp_id
    from_index_to_BaPCodindex(v_id, var_bcid)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    # Register the variable
    c_register_variable(m.ptr, name, Cint(column_id - 1), sp_type, sp_bcid, var_bcid)
  end
  c_init_variables(m.ptr, l, u, c)
end

# Entry method for registration without variable names and indexes
function c_register_vars(m::BcModel, A, l, u, c)
  nb_rows, nb_cols = size(A)
  empty = [0, -1, -1, -1, -1, -1, -1, -1]
  for col_id in 1:nb_cols
    c_register_variable(m.ptr, :var, Cint(col_id - 1), :MIP, empty, empty)
  end
  c_init_variables(m.ptr, l, u, c)
end

# Registering a variable in BaPCod
function c_register_variable(mptr::Ptr{Cvoid}, name::Symbol, column_id::Cint,
                             sp_type::Symbol, sp_mid::Vector{Int}, var_mid::Vector{Int})
  sp_mid2 = [Cint(s) for s in sp_mid]
  var_mid2 = [Cint(s) for s in var_mid]
  c_register_variable(mptr, name, column_id, sp_type, sp_mid2, var_mid2)
end

function c_register_variable(mptr::Ptr{Cvoid}, name::Symbol, column_id::Cint,
                             sp_type::Symbol, sp_mid::Vector{Cint}, var_mid::Vector{Cint})
  if sp_type != :ALL
    status = register_var!(mptr, name, column_id, sptype_to_int(sp_type), sp_mid, var_mid)
  else
    status = register_generic_var!(mptr, name, column_id)
  end
  (status != 1) && error("Cannot init variable $name at column $column_id.")
end

# Init variables in BaPCod (store bounds and costs)
function c_init_variables(mptr::Ptr{Cvoid}, l, u, c)
  init_vars!(mptr, l, u, c)
end

function c_set_sp_multiplicities(m::BcModel, sp_mult)
  sp_bcid = Array{Cint}(undef, 8)
  for (sp_id, sp_type, mult_lb, mult_ub) in sp_mult
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    status = sub_problem_mult!(m.ptr, Cint(mult_lb), Cint(mult_ub), sptype_to_int(sp_type), sp_bcid)
    (status != 1) && error("Cannot set multiplicity on the subproblem with the index $sp_id. Make sure it exists.")
  end
end

function c_set_sp_priorities(m::BcModel, sp_prio)
  sp_bcid = Array{Cint}(undef, 8)
  for (sp_id, priority) in sp_prio
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    status = set_subproblem_priority!(m.ptr, sp_bcid, priority)
    (status != 1) && error("Cannot set priority on the subproblem with the index $sp_id. Make sure it exists.")
  end
end

# Not called by (Block)JuMP
function c_get_var_lb(model::BcModel)
  nbvars = nb_vars(model)
  lb = Array{Cdouble}(undef, nbvars)
  status = get_var_lb!(model.ptr, lb, nbvars)
  (status != 1) && error("Cannot get variables lower bounds.")
  lb
end

function c_set_var_lb!(model::BcModel, l)
  error("Function bc_vars.jl > c_set_var_lb! can not be used because BaPCod_framework may not end the resolution properly (all containers cleaned...). You should rather solve a new JuMP.Model")
  size = length(l)
  status = set_var_lb!(model.ptr, l, size)
  (status != 1) && error("Cannot set variables lower bounds.")
end

# Not called by (Block)JuMP
function c_get_var_ub(model::BcModel)
  nbvars = nb_vars(model)
  ub = Array{Cdouble}(undef, nbvars)
  status = get_var_ub!(model.ptr, ub, nbvars)
  (status != 1) && error("Cannot get variables lower bounds.")
  ub
end

function c_set_var_ub!(model::BcModel, u)
  size = length(u)
  status = set_var_ub!(model.ptr, u, size)
  (status != 1) && error("Cannot set variables upper bounds.")
end

# Not called by (Block)JuMP
function c_get_var_costs(model::BcModel)
  nbvars = nb_vars(model)
  c = Array{Cdouble}(undef, nbvars)
  status = get_var_costs!(model.ptr, c, nbvars)
  (status != 1) && error("Cannot get variables costs.")
  c
end

function c_set_var_costs!(model::BcModel, c)
  size = length(c)
  status = set_var_costs!(model.ptr, c, size)
  (status != 1) && error("Cannot set variables costs.")
end

const rev_var_type_map = Dict(
  :Cont => 'C',
  :Bin => 'B',
  :Int => 'I',
)

# Not called by (Block)JuMP
function c_get_var_type(model::BcModel)
  nbvars = nb_vars(model)
  t = Array{Cchar}(undef, nbvars)
  status = get_var_type!(model.ptr, t, nbvars)
  (status != 1) && error("Cannot get variables types.")
  t
end

function c_set_var_type!(model::BcModel, typ::Vector{Symbol})
  bctyp = join(map(x -> rev_var_type_map[x], typ))
  size = length(typ)
  status = set_var_type!(model.ptr, bctyp, size)
  (status != 1) && error("Cannot set variables types.")
end

function c_vars_branching_priorities(model::BcModel, p::Dict{Tuple{Symbol, Tuple, Symbol}, Cdouble})
  sp_bcid = Array{Cint}(undef, 8)
  for ((varname, (sp_type, sp_id), where_), priority) in p
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    if where_ == :MASTER
      status = set_var_priority_in_master!(model.ptr, varname, sp_bctype, sp_bcid, priority)
    else
      status = set_var_priority_in_sp!(model.ptr, varname, sp_bctype, sp_bcid, priority)
    end
    (status != 1) && error("Cannot set branching priority on variables named $varname.")
  end
end

function c_branching_rules(model::BcModel, bc_branchings::Dict{Symbol, Any})
  c_model = model.ptr
  for branching in bc_branchings
    if branching[1] == :aggrsubprob_variables
      for (varname, args) in branching[2]
        # Default values of parameters
        highest_priority = 0.5
        priority = 1
        number_of_ignored_indices = 0
        preprocessing = false
        # User parameters values
        for (key, val) in args
          if key == :highest_priority
            highest_priority = val
          elseif key == :priority
            priority = val
          elseif key == :number_of_ignored_indices
            number_of_ignored_indices = val
          elseif key == :preprocessing
            preprocessing = val
          else
            error("Unrecognized parameter $key for branching rule $(branching[1]).")
          end
        end
        #info("branching rule $(branching[1]) on variable $varname : highest_priority = $highest_priority -- priority = $priority -- number_of_ignored_indices = $number_of_ignored_indices -- preprocessing = $preprocessing.")
        wbcm_add_aggrsubprob_var_branching(c_model, varname, Cdouble(highest_priority), Cdouble(priority), Cint(number_of_ignored_indices), preprocessing)
      end
    else
      error("Unrecognized BaPCod branching rule $(branching[1]).")
    end
  end
end

function c_branching_exp(model::BcModel, exp_dict::Dict)
   expbcid = Array{Cint}(undef, 8)
   for (name, exp) in exp_dict
    p = exp[1][4]
    arrayid = wbcm_register_branching_expression(model.ptr, name, p) 
    for (index, coeff, var, p) in exp
            from_index_to_BaPCodindex(index, expbcid)
            varsid = Cint.(var .- 1)
            wbcm_add_branching_expression(model.ptr, arrayid, expbcid, varsid, coeff, Cint(length(var)))
        end
   end
   return
end
