abstract type BcCallbackData <: MathProgCallbackData end

mutable struct BcSepCallbackData <: BcCallbackData
  solptr::Ptr{Cvoid}
  model::BcMathProgModel
  cutlist::Ptr{Cvoid}
  idcb::Int
  state::Symbol
end

function mastercallback(solptr::Ptr{Cvoid}, userdata::Ptr{Cvoid}, cutlist::Ptr{Cvoid}, idfunctor::Cint)
  model = unsafe_pointer_to_objref(userdata)::BcMathProgModel
  # Lazy callback
  if model.lazycb != nothing
    bccbd = BcSepCallbackData(solptr, model, cutlist, idfunctor, :Unknown)
    status = model.lazycb(bccbd)
    if status === :Exit
      println("Stop the solver. Will quit().")
      quit()
    end
  end
  # User callback
  if model.usercb != nothing
    bccbd = BcSepCallbackData(solptr, model, cutlist, idfunctor, :Unknown)
    status = model.usercb(bccbd)
    if status === :Exit
      println("Stop the solver. Will quit().")
      quit()
    end
  end
  nothing
end

function register_lazycb!(m::BcMathProgModel)
  c_model = m.inner.ptr
  lazycb_c_fct = @cfunction(mastercallback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint))
  status = wbcs_init_sep_routine!(c_model, lazycb_c_fct, 'C', m)
  (status != 1) && error("Cannot initialize the lazy cut callback.")
end

function register_usercb!(m::BcMathProgModel)
  c_model = m.inner.ptr
  usercb_c_fct = @cfunction(mastercallback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint))
  status = wbcs_init_sep_routine!(c_model, usercb_c_fct, 'F', m)
  (status != 1) && error("Cannot initialize the user cut callback.")
end

getsolutioncost(d::BcSepCallbackData) = getsolutioncost(d.solptr)
getsolutiontruecost(d::BcSepCallbackData) = getsolutiontruecost(s.solptr)

# Adds cut callback f to the model. Function f takes as argument only a
# MathProgCallbackData object.
function setcutcallback!(m::BcMathProgModel, f)
  m.usercb = f
end

# Adds lazy cut callback f to the model
function setlazycallback!(m::BcMathProgModel, f)
  m.lazycb = f
end

# Adds lazy constraint callback f to the model.
# Function f takes as argument only a MathProgCallbackData object.
function setinfocallback!(m::BcMathProgModel, f)
  error("setinfocallback! not yet implemented")
end

# Grabs current best integer-feasible solution to the model.
# The optional second argument specifies an output vector.
function cbgetmipsolution(d::BcCallbackData)
  @assert d.state == :MIPSol
  sol::Array{Cdouble}(nb_vars(d.model))
  getaggregatedsolution(d.model.inner, d.solptr, sol)
  return sol
end

function cbgetmipsolution(d::BcCallbackData, sol::Vector{Cdouble})
  @assert d.state == :MIPSol
  getaggregatedsolution(d.model.inner, d.solptr, sol)
  return nothing
end

function cbgetlpsolution(d::BcCallbackData)
  @assert d.state == :MIPNode || d.state == :MIPBranch
  sol::Array{Cdouble}(nb_vars(d.model))
  getaggregatedsolution(d.model.inner, d.solptr, sol)
  return sol
end

function cbgetlpsolution(d::BcCallbackData, sol::Vector{Cdouble})
  @assert d.state == :MIPNode || d.state == :MIPBranch
  getaggregatedsolution(d.model.inner, d.solptr, sol)
  return nothing
end

# ##############################################################################
# Returns current location in solve process:
# - :MIPNode if at node in branch-and-cut tree,
# - :MIPSol at an integer-feasible solution,
# - and :Intermediate otherwise.
# ##############################################################################
# - MIPNode: when we are at a node in the branch-and-cut tree. This is generally
# used for access to the solution of a relaxation (via cbgetlpsolution()), or
# for adding cuts (via cbaddcut!() or cbaddcutlocal!()).
# - MIPSol: when we have found a new MIP incumbent (an integer-feasible solution).
# This is generally used for keeping track of the intermediate solutions
# generated (via cbgetmipsolution()) inside the branch-and-cut tree.
# - Intermediate: when we are still in the process of MIP or during iterations
# of a continuous solver. For MIPs, this is generally be used for keeping track
# of pessimistic and optimistic bounds (via cbgetobj() and cbgetbestbound()
# respectively), or the number of explored nodes (via cbgetexplorednodes()) in
# the branch-and-cut tree.
# ##############################################################################
# BaPCod knows when he must add lazy & user cuts
function cbgetstate(d::BcCallbackData)
  cbtype = c_get_sep_cb_type!(d.model.inner.ptr, Cint(d.idcb))

  if cbtype == 67 # C as Core = Lazy
    d.state = :MIPSol
  elseif cbtype == 70 # F as Facultative = User
    d.state = :MIPNode
  else
    d.state = :Intermediate
  end
  return d.state
end

#  Adds cut to model. The coefficient values are represented sparsely, with
# (one-indexed) indices in varidx and values in varcoef. The constraint sense
# sense is a character taking value <, >, or =, and the right-hand side value
# is rhs.
function cbaddcut!(d::MathProgCallbackData, varidx, varcoef, sense, rhs)
  c_model = d.model.inner.ptr
  nb_vars = length(varcoef)
  status = wbcs_add_sep_cut!(c_model, d.idcb, d.cutlist, varcoef, varidx,
                             nb_vars, sense, rhs)
  (status != 1) && error("Cannot add the user cut to the model.")
end

# Adds local cut to model. It works as cbaddcut! but the cut is local in the
# sense that it only applies to the current node and the subtree rooted at this
# node.
cbaddcutlocal!(d::MathProgCallbackData, varidx, varcoef, sense, rhs) =
  error("Local cuts are not supported by BaPCod.")

# Adds lazy constraint to model. The coefficient values are represented sparsely,
# with (one-indexed) indices in varidx and values in varcoef. The constraint
# sense sense is a character taking value <, >, or =, and the right-hand side
# value is rhs.
function cbaddlazy!(d::MathProgCallbackData, varidx, varcoef, sense, rhs)
  c_model = d.model.inner.ptr
  nb_vars = length(varcoef)
  status = wbcs_add_sep_cut!(c_model, Cint(d.idcb), d.cutlist, varcoef, varidx,
                             Cint(nb_vars), sense, rhs)
  (status != 1) && error("Cannot add the lazy cut to the model.")
end

# Adds local lazy constraint to model. It works as cbaddlazy! but the lazy
# constraint is local in the sense that it only applies to the current node and
# the subtree rooted at this node.
cbaddlazylocal!(d::MathProgCallbackData, varidx, varcoef, sense, rhs) =
  error("Local cuts are not supported by BaPCod.")


# Submit a (possibly partially defined) heuristic solution for the model. Should
# reset the solution stored in d to the original state at the start of callback.
function cbaddsolution!(d::MathProgCallbackData)
  error("cbaddsolution not yet implemented")
end

# Sets the value of a variable with (one-based) index varidx to value in the
# current partial solution being constructed by a user heuristic.
function cbsetsolutionvalue!(d::MathProgCallbackData, varidx, value)
  error("cbsetsolutionvalue not yet implemented")
end
