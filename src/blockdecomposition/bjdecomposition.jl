function add_Dantzig_Wolfe_decomposition(m::JuMP.Model, f::Function)
  m.ext[:block_decomposition].DantzigWolfe_decomposition_fct = f
end


function add_Dantzig_Wolfe_decomposition_on_variables(m::JuMP.Model, f::Function)
  m.ext[:block_decomposition].DantzigWolfe_decomposition_on_vars_fct = f
end

function use_DantzigWolfe(m::JuMP.Model)
  return m.ext[:block_decomposition].DantzigWolfe_decomposition_fct != nothing
end

"""
    add_Benders_decomposition(model::JuMP.Model, B_decomp::Function)

assigns the decomposition function `B_decomp` to the `model`.
"""

function add_Benders_decomposition(m::JuMP.Model, f::Function)
  m.ext[:block_decomposition].Benders_decomposition_fct = f
end

function add_Benders_decomposition_on_constraints(m::JuMP.Model, f::Function)
  error("Not supported")
  m.ext[:block_decomposition].Benders_decomposition_on_cstrs_fct = f
end

function use_Benders_decomposition(m::JuMP.Model)
  return m.ext[:block_decomposition].Benders_decomposition_fct != nothing
end

function create_cstrs_vars_decomposition_list!(m::JuMP.Model)
  A = prepConstrMatrix(m)
  m.ext[:cstrs_decomposition_list] = create_cstrs_decomposition_list(m, sparse(A'))
  m.ext[:vars_decomposition_list] = create_vars_decomposition_list(m, A)
end

function ismip(m::JuMP.Model)
  DW_dec_f = m.ext[:block_decomposition].DantzigWolfe_decomposition_fct
  B_dec_f = m.ext[:block_decomposition].Benders_decomposition_fct
  ret = isempty(m.ext[:generic_vars]) && isempty(m.ext[:generic_cstrs]) && DW_dec_f == nothing && B_dec_f == nothing
  return ret
end

function DW_decomposition(DW_dec_f, cstr_name, cstr_id)
  (sp_type, sp_id) = (nothing, nothing)
  if applicable(DW_dec_f, cstr_name, cstr_id)
    (sp_type, sp_id) = DW_dec_f(cstr_name, cstr_id)
    if (sp_type != :DW_MASTER && sp_type != :DW_SP) || (!isa(sp_id, Tuple) && !isa(sp_id, Integer))
      info = """ The Dantzig-Wolfe decomposition function must return a tuple (sp_type::Symbol, sp_id::Union{Integer,Tuple})
                 - sp_type is the type of the problem to which the constraint is assigned
                   it must be :DW_MASTER for the master or :DW_SP for the subproblem
                 - sp_id is a the index of the subproblem.
             """
      bjerror(info, "Dantzig-Wolfe decomposition function does not return a correct tuple.")
    end
  else
    info = """ The Dantzig-Wolfe decomposition function must take the two following arguments :
               - cstr_name::Symbol, the name of the constraint
               - cstr_id::Union{Integer, Tuple}, the index of the constraint.
           """
    errmsg = "Cannot applicate the function with the arguments $cstr_name::$(typeof(cstr_name)) and $cstr_id::$(typeof(cstr_id))"
    bjerror(info, errmsg)
  end
  if isa(sp_id, Integer)
    sp_id = (sp_id,)
  end
  (sp_type, sp_id)
end

function B_decomposition(B_dec_f, var_name, var_id)
  (sp_type, sp_id) = (nothing, nothing)
  if applicable(B_dec_f, var_name, var_id)
    (sp_type, sp_id) = B_dec_f(var_name, var_id)
    if (sp_type != :B_MASTER && sp_type != :B_SP) || (!isa(sp_id, Tuple) && !isa(sp_id, Integer))
      info = """ The Benders decomposition function must return a tuple (sp_type::Symbol, sp_id::Union{Integer, Tuple})
                 - sp_type is the type of the problem to which the constraint is assigned
                   it must be :B_MASTER for the master or :B_SP for the subproblem
                 - sp_id is a the index of the subproblem.
             """
      errmsg = "Benders decomposition function does not return a correct tuple ($sp_type::$(typeof(sp_type)), $sp_id::$(typeof(sp_id)))."
      bjerror(info, errmsg)
    end
  else
    info = """ The Benders decomposition function must take the two following arguments :
               - var_name::Symbol, the name of the variable
               - var_id::Union{Integer, Tuple}, the index of the variable.
           """
    errmsg = "Cannot applicate the function with the arguments $var_name::$(typeof(var_name)) and $var_id::$(typeof(var_id))"
    bjerror(info, errmsg)
  end
  if isa(sp_id, Integer)
    sp_id = (sp_id,)
  end
  (sp_type, sp_id)
end

function create_cstrs_decomposition_list(m::JuMP.Model, A)
  sp_id = 0
  sp_type = :MIP
  rows = rowvals(A)
  DW_dec_f = m.ext[:block_decomposition].DantzigWolfe_decomposition_fct
  B_dec_f = m.ext[:block_decomposition].Benders_decomposition_fct
  cstrs_list = Array{Tuple}(undef, size(m.ext[:varcstr_report].cstrs_report))
  mip = ismip(m)

  for (row_id, (name, cstr_id)) in enumerate(m.ext[:varcstr_report].cstrs_report)
    if !mip
      sp_type = :DW_MASTER
      if is_genericcstr(m, row_id) #Is it a generic constraint ?
        sp_type = :ALL
      elseif DW_dec_f != nothing # Dantzig-Wolfe decomposition
        if name == :anonymous_cstr
          #@warn("Anonymous constraint assigned to master problem during decomposition.")
          sp_type = :DW_MASTER
        else
          (sp_type, sp_id) = DW_decomposition(DW_dec_f, name, cstr_id)
        end
      elseif B_dec_f != nothing # Benders decomposition
        sp_type = :B_MASTER
        nb_vars = 0
        for i in nzrange(A, row_id)
          (var_name, var_id) = m.ext[:varcstr_report].vars_report[rows[i]]
          (var_sp_type, var_sp_id) = B_decomposition(B_dec_f, var_name, var_id)
          if nb_vars == 0 || (nb_vars > 0 && var_sp_type != :B_MASTER)
            # Check if the constraint is in the same subproblem
            if nb_vars > 0 && sp_type != :B_MASTER && var_sp_id != sp_id
              bjerror("A single constraint cannot belongs to two different subproblems. (constraint = $name & id = $cstr_id)")
            end
            sp_type = var_sp_type
            sp_id = var_sp_id
          end
          nb_vars += 1
        end
      end
    end
    cstrs_list[row_id] = (name, cstr_id, sp_type, sp_id)
    list_sp!(m, sp_type, sp_id)
  end
  cstrs_list
end

function is_genericvar(m::JuMP.Model, varname::Symbol)
  return haskey(m.ext[:generic_vars], varname)
end

function is_genericcstr(m::JuMP.Model, cstrid::Int)
  return haskey(m.ext[:generic_cstrs], cstrid)
end

function create_vars_decomposition_list(m::JuMP.Model, A)
  sp_id = 0
  sp_type = :MIP
  rows = rowvals(A)
  DW_dec_f = m.ext[:block_decomposition].DantzigWolfe_decomposition_fct
  DW_dec_on_vars_f = m.ext[:block_decomposition].DantzigWolfe_decomposition_on_vars_fct
  B_dec_f = m.ext[:block_decomposition].Benders_decomposition_fct
  vars_list = Array{Tuple}(undef, size(m.ext[:varcstr_report].vars_report))
  mip = ismip(m)

  for (column_id, (name, var_id)) in enumerate(m.ext[:varcstr_report].vars_report)
    if !mip
      sp_type = :DW_MASTER
      if is_genericvar(m, name) #Is it a generic variable ?
        sp_type = :ALL
      elseif B_dec_f != nothing # benders decomposition
        (sp_type, sp_id) = B_decomposition(B_dec_f, name, var_id)
      elseif DW_dec_f != nothing # Dantzig-Wolfe decomposition
        if DW_dec_on_vars_f != nothing
          (sp_type, sp_id) = DW_decomposition(DW_dec_on_vars_f, name, var_id)
        else
          sp_type = :DW_MASTER # If the variable is not used in cstrs => go to master
          nb_cstrs = 0
          for j in nzrange(A, column_id)
            (cstr_name, cstr_id) = m.ext[:varcstr_report].cstrs_report[rows[j]]
            (cstr_sp_type, cstr_sp_id) = DW_decomposition(DW_dec_f, cstr_name, cstr_id)
            if (nb_cstrs == 0) || (nb_cstrs > 0 && cstr_sp_type != :DW_MASTER)
              # Check if the variable is in the same subproblem (except MASTER)
              if nb_cstrs > 0 && sp_type != :DW_MASTER && cstr_sp_id != sp_id
                bjerror("A single variable cannot belong to two different subproblems. (variable = $name & id = $var_id).")
              end
              sp_type = cstr_sp_type
              sp_id = cstr_sp_id
            end
            nb_cstrs += 1
          end
        end
      end
    end
    vars_list[column_id] = (name, var_id, sp_type, sp_id)
    list_sp!(m, sp_type, sp_id)
  end
  vars_list
end

function list_sp!(m::JuMP.Model, sp_type, sp_id)
  if sp_type == :DW_SP
    if !haskey(m.ext[:sp_list_dw], sp_id)
      m.ext[:sp_list_dw][sp_id] = 1
    end
  elseif sp_type == :B_SP
    if !haskey(m.ext[:sp_list_b], sp_id)
      m.ext[:sp_list_b][sp_id] = 1
    end
  end
end

function create_sp_mult_tab!(m::JuMP.Model)
  if m.ext[:sp_mult_fct] != nothing
    m.ext[:sp_mult_tab] = Array{Tuple}(undef, 0)
    fill_sp_mult_tab!(m, :sp_list_dw, :DW_SP)
    fill_sp_mult_tab!(m, :sp_list_b, :B_SP)
  end
end

function fill_sp_mult_tab!(m::JuMP.Model, tab::Symbol, sp_type::Symbol)
  if length(m.ext[tab]) > 0
    mult_fct = m.ext[:sp_mult_fct]
    fkspid = first(m.ext[tab])[1] # Get the key of the first entry
    # Check is the function is applicable only on the first subproblem
    if !applicable(mult_fct, fkspid, sp_type)
      info = """First argument = $fkspid (type = $(typeof(fkspid)))
                Second argument = $sp_type (type = $(typeof(fkspid))) """
      bjerror(info, "The function defining multiplicities of subproblems is not applicable.")
    end
    for sp_id in keys(m.ext[tab])
      (lb, ub) = mult_fct(sp_id, sp_type)
      if !isa(lb, Integer) || !isa(ub, Integer)
        bjerror("You use lb = $lb and ub = $ub", "Multiplicity must be an integer.")
      end
      if ub < lb
        bjerror("Multiplicity (ub < lb) : ($ub < $lb)")
      end
      push!(m.ext[:sp_mult_tab], (sp_id, sp_type, lb, ub))
    end
  end
end

function create_sp_prio_tab!(m::JuMP.Model)
  if m.ext[:sp_prio_fct] != nothing
    m.ext[:sp_prio_tab] = Array{Tuple}(0)
    fill_sp_prio_tab!(m, :sp_list_dw, :DW_SP)
    fill_sp_prio_tab!(m, :sp_list_b, :B_SP)
  end
end

function fill_sp_prio_tab!(m::JuMP.Model, tab::Symbol, sp_type::Symbol)
  if length(m.ext[tab]) > 0
    priority_fct = m.ext[:sp_prio_fct]
    fkspid = first(m.ext[tab])[1] # Get the key of the first entry
    # Check is the function is applicable only on the first subproblem
    if !applicable(priority_fct, fkspid, sp_type)
      info = """First argument = $fkspid (type = $(typeof(fkspid)))
                Second argument = $sp_type (type = $(typeof(fkspid))) """
      bjerror(info, "The function defining priorities of subproblems is not applicable.")
    end
    for sp_id in keys(m.ext[tab])
      priority = priority_fct(sp_id, sp_type)
      if !isa(priority, Integer) || priority < 0 || priority >= Inf
        bjerror("You use $priority.", "Priority must be a positive integer.")
      end
      push!(m.ext[:sp_prio_tab], (sp_id, sp_type, priority))
    end
  end
end
