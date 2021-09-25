mutable struct BcGenVarCallbackData <: BcCallbackData
  memberships::Ptr{Cvoid}  # List of (cstr, coeff). The variable will be added in
                         # each constraint cstr with a coefficient coeff
  objcoeff
  varindex::Tuple
  variable::JuMP.Variable
end

mutable struct BcGenCstrCallbackData <: BcCallbackData
  constraint::JuMP.ConstraintRef
  bcconstraint::Ptr{Cvoid}
  cstrindex::Tuple
end

getcstrindex(cb::BcGenCstrCallbackData) = cb.cstrindex
getvarindex(cb::BcGenVarCallbackData) = cb.varindex

function genvarcallback(index::Ptr{Cint}, var::Ptr{Cvoid}, cstrsmembership::Ptr{Cvoid}, obj::Ptr{Cdouble})
  variable = unsafe_pointer_to_objref(var)::JuMP.Variable # get the JuMP.Variable
  varindex = from_BaPCodindex_to_index(unsafe_wrap(Array, index, 8; own = false)) # get the variable index
  objcoeff = unsafe_wrap(Array, obj, 1; own = false)
  cbdata = BcGenVarCallbackData(cstrsmembership, objcoeff, varindex, variable)
  # Callback
  internalmodel(variable.m).solver.genvarcallbacks[Symbol(getname(variable))][2](cbdata)
  return
end

# Callbacks for generic variables
function register_genvarcb!(m::BcMathProgModel)
  for genvarcb in m.solver.genvarcallbacks
    name = genvarcb[1]
    jump_var, fct = genvarcb[2]
    c_fct = @cfunction(genvarcallback, Cvoid, (Ptr{Cint}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
    status = wbcm_attach_var_functor!(m.inner.ptr, name, c_fct, jump_var)
    (status != 1) && error("Cannot attach the generic variable callback to the variable $var.")
  end
  return
end

function gencstrcallback(index::Ptr{Cint}, cstr::Ptr{Cvoid}, bccstr::Ptr{Cvoid})
  constraint = unsafe_pointer_to_objref(cstr)::JuMP.ConstraintRef
  cstrindex = from_BaPCodindex_to_index(unsafe_wrap(Array, index, 8; own = false))
  cbdata = BcGenCstrCallbackData(constraint, bccstr, cstrindex)
  internalmodel(constraint.m).solver.gencstrcallbacks[constraint.idx][3](cbdata)
  return
end

function register_gencstrcb!(m::BcMathProgModel)
  for gencstrcb in m.solver.gencstrcallbacks
    name = gencstrcb[1]
    jump_cstr, name, fct = gencstrcb[2]
    c_fct = @cfunction(gencstrcallback, Cvoid, (Ptr{Cint}, Ptr{Cvoid}, Ptr{Cvoid}))
    status = wbcm_attach_cstr_functor!(m.inner.ptr, Symbol(name), c_fct, jump_cstr)
    (status != 1) && error("Cannot attach the functor to the constraint")
  end
  return
end

function adddynamicvariable!(m::BcMathProgModel, name::String, id::Tuple, col_id::Int, sp_type::Symbol, sp_id::Tuple)
  c_model = m.inner.ptr
  col_id = Cint(col_id - 1)
  var_bcid = from_index_to_BaPCodindex(id)
  sp_bctype = Cint(sptype_to_int(sp_type))
  sp_bcid = from_index_to_BaPCodindex(sp_id)
  status = wbcm_add_dyn_var!(c_model, col_id, Symbol(name), var_bcid, sp_bctype, sp_bcid)
  (status != 1) && error("Cannot generate the variable with id $id.")
  return
end

function adddynamiccut!(m::BcMathProgModel, name::String, id::Tuple, row_id::Int, cuts_list::Ptr{Cvoid})
  c_model = m.inner.ptr
  row_id = Cint(row_id - 1)
  cstr_bcid = from_index_to_BaPCodindex(id)
  status = wbcm_add_dyn_cstr!(c_model, row_id, Symbol(name), cstr_bcid)
  (status != 1) && error("Cannot generate the constraint with id $id.")
  cstr = wbcm_get_cstr(c_model, row_id)
  (cstr == C_NULL) && error("Cannot get the generated constraint")
  wbcs_add_cut!(cuts_list, cstr)
  return
end

function adddynamicconstraint!(m::BcMathProgModel, name::String, id::Tuple, row::Int, cb::BcOracleSolverData)
  c_model = m.inner.ptr
  row_id = Cint(row - 1)
  cstr_bcid = from_index_to_BaPCodindex(id)
  status = wbcm_add_dyn_cstr!(c_model, row_id, Symbol(name), cstr_bcid)
  (status != 1) && error("Cannot generate the constraint with id $id.")
  cstr = wbcm_get_cstr(c_model, row_id)
  (cstr == C_NULL) && error("Cannot get the generated constraint")
  status = wbcm_set_constr_bc_type!(cstr, 'S')
  (status != 1) && error("Cannot set the type of the constraint.")
  wbcs_add_to_oracle_dual_sol!(c_model, cb.solptr, cstr)
  return
end

function addtermtodynamiccstr!(cb::BcGenCstrCallbackData, col_id::Cint, coeff)
  c_model = internalmodel(cb.constraint.m).inner.ptr
  col_id = Cint(col_id - 1)
  status = wbcm_add_cstr_term!(c_model, cb.bcconstraint, col_id, Cdouble(coeff))
  (status != 1) && error("Cannot add a term to the constraint.")
  return
end

Base.cconvert(::Type{Ptr{Cint}}, array::Vector{Int}) = [Cint(i) for i in array]

function addtermstodynamiccstr!(cb::BcGenCstrCallbackData, cols_id::Vector{Int}, coeff)
  c_model = internalmodel(cb.constraint.m).inner.ptr
  status = wbcm_add_cstr_terms!(c_model, cb.bcconstraint, Base.cconvert(Ptr{Cint}, cols_id), Cdouble(coeff))
  (status != 1) && error("Cannot add terms to the constraint.")
  return
end

function addentrytodynamiccol!(cb::BcGenVarCallbackData, row_id::Integer, coeff)
  c_model = internalmodel(cb.variable.m).inner.ptr
  row_id = Cint(row_id - 1)
  status = wbcm_add_membership_to_cstr!(c_model, cb.memberships, row_id, Cdouble(coeff))
  (status != 1) && error("Cannot add en entry to the column.")
  return
end

function setobjvaluetodynamiccol!(cb, value)
  cb.objcoeff[1] = value
  return
end

function addrhstodynamiccstr!(cb::BcGenCstrCallbackData, sense::Char, rhs)
  wbcm_add_cstr_rhs!(cb.bcconstraint, sense, Cdouble(rhs))
  return
end
