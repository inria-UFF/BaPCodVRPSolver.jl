mutable struct BcSepCutCallbackData <: BcCallbackData # TODO rename 
  solptr::Ptr{Cvoid}
  model::BcMathProgModel
  cutlist::Ptr{Cvoid}
  idcb::Int
  state::Symbol
  sol::Vector{Cdouble}
end

function master_corecuts_callback(solptr::Ptr{Cvoid}, userdata::Ptr{Cvoid}, cutlist::Ptr{Cvoid}, idfunctor::Cint)::Cvoid
  model = unsafe_pointer_to_objref(userdata)::BcMathProgModel
  bccbd = BcSepCutCallbackData(solptr, model, cutlist, idfunctor, :Unknown, Vector{Cdouble}())
  cbgetsolution(bccbd)
  model.solver.corecuts_cb_func(bccbd)
  nothing
end

function master_facultativecuts_callback(solptr::Ptr{Cvoid}, userdata::Ptr{Cvoid}, cutlist::Ptr{Cvoid}, idfunctor::Cint)::Cvoid
  model = unsafe_pointer_to_objref(userdata)::BcMathProgModel
  bccbd = BcSepCutCallbackData(solptr, model, cutlist, idfunctor, :Unknown, Vector{Cdouble}())
  cbgetsolution(bccbd)
  model.solver.facultativecuts_cb_func(bccbd)
  nothing
end

function register_core_cb!(m::BcMathProgModel)
  c_model = m.inner.ptr
  lazycb_c_fct = cfunction(master_corecuts_callback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint))
  status = wbcs_init_sep_routine!(c_model, lazycb_c_fct, 'C', m)
  (status != 1) && error("Cannot initialize the lazy cut callback.")
end

function register_facultative_cb!(m::BcMathProgModel)
  c_model = m.inner.ptr
  usercb_c_fct = cfunction(master_facultativecuts_callback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint))
  status = wbcs_init_sep_routine!(c_model, usercb_c_fct, 'F', m)
  (status != 1) && error("Cannot initialize the user cut callback.")
end

function cbgetsolution(d::BcCallbackData)
  sol = Array{Cdouble}(nb_vars(d.model))
  getaggregatedsolution(d.model.inner, d.solptr, sol)
  d.sol = copy(sol)
end

function getvalvariable(d::BcCallbackData, varidx)
  return d.sol[varidx]
end
