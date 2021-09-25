mutable struct BcOracleSolverData <: OracleSolverData
  createnewsol::Bool
  solptr::Ptr{Cvoid}
  model::BcMathProgModel
  sp_id::Tuple
  sp_type::Symbol
end

function register_oracles!(m::BcMathProgModel)
  for (oracle_id, (sp_id, sp_type, oracle_fct)) in enumerate(m.solver.oracles)
    sp_bcid = Vector{Cint}(8)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    oraclesolver_fct = cfunction(oraclecallback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Cint))
    status = c_initOracle!(m.inner.ptr, oraclesolver_fct, sptype_to_int(sp_type), sp_bcid, m, oracle_id)
    (status != 1) && error("Cannot assign oracle to subproblem $sp_id of type $sp_type.")
  end
end

function oraclecallback(solptr::Ptr{Cvoid}, userdata::Ptr{Cvoid}, idoracle::Cint)
  model = unsafe_pointer_to_objref(userdata)::BcMathProgModel
  oracle = model.solver.oracles[idoracle]
  bcosd = BcOracleSolverData(false, solptr, model, oracle[1], oracle[2])
  oracle[3](bcosd) # call of the function defined by the user
  nothing
end

set_oracles!(s::BaPCodSolver, oracles) = (s.oracles = oracles)

function set_oraclesolution_newsolution(o::BcOracleSolverData)
  o.createnewsol = true
end

function set_oraclesolution_solution(o::BcOracleSolverData, x::JuMP.Variable, v::Real)
  if o.createnewsol
    status = c_newOracleSol!(o.solptr)
    (status != 1) && error("Cannot attach a new solution to the oracle solution.")
    o.createnewsol = false
  end
  s = c_addToOracleSol(o.model.inner.ptr,o.solptr,Cint(x.col-1),Cdouble(v))
  (s != 1) && error("Cannot add new solution ($x = $v) to the oracle solution.")
end

function set_oraclesolution_solution(o::BcOracleSolverData, x::BlockDecompositionExtras.DynVariable, v::Real)
  if o.createnewsol
    status = c_newOracleSol!(o.solptr)
    (status != 1) && error("Cannot attach a new solution to the oracle solution.")
    o.createnewsol = false
  end
  s = c_addToOracleSol(o.model.inner.ptr,o.solptr,Cint(x.col-1),Cdouble(v))
  (s != 1) && error("Cannot add new solution ($x = $v) to the oracle solution.")
end

function set_oraclesolution_objval(o::BcOracleSolverData, objval::Real)
  s = c_setObjValOfOracleSol(o.solptr, Cdouble(objval))
  (s != 1) && error("Cannot set the objective value of the oracle solution.")
end

get_oracle_phaseofstageapproach(o::BcOracleSolverData) =
  c_getOraclePhaseOfStage(o.solptr)

set_solution_to_homogeneous_syst(o::BcOracleSolverData) =
  c_setSolToHomogeneousSyst(o.solptr)

set_solution_to_homogeneous_syst_with_eq(o::BcOracleSolverData) =
  c_setSolToHomogeneousSystWithEq(o.solptr)
