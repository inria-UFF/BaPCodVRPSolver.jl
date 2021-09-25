import MathProgBase.SolverInterface:
  LinearQuadraticModel,
  loadproblem!,
  writeproblem,
  getvarLB,
  setvarLB!,
  getvarUB,
  setvarUB!,
  getconstrLB,
  setconstrLB!,
  getconstrUB,
  setconstrUB!,
  getobj,
  setobj!,
  getconstrmatrix,
  addvar!,
  addconstr!,
  numlinconstr,
  getconstrsolution,
  getreducedcosts,
  getconstrduals,
  getsimplexiter,
  setwarmstart!,
  getinfeasibilityray,
  getbasis,
  getunboundedray,
  getbarrieriter,
  setcutcallback!,
  setlazycallback!,
  setinfocallback!,
  cbgetmipsolution,
  cbgetlpsolution,
  cbgetstate,
  cbaddcut!,
  cbaddcutlocal!,
  cbaddlazy!,
  cbaddlazylocal!,
  cbaddsolution!,
  cbsetsolutionvalue!,
  getsolution,
  getobjval,
  optimize!,
  status,
  getobjbound,
  getobjgap,
  getrawsolver,
  getsolvetime,
  getnodecount,
  setsense!,
  getsense,
  numvar,
  numconstr,
  freemodel!,
  setvartype!,
  getvartype,
  setparameters!

import ...BlockDecomposition.BlockSolverInterface:
  set_constrs_decomposition!,
  set_vars_decomposition!,
  set_sp_mult!,
  set_var_branching_prio!,
  set_sp_ids!,
  set_sp_prio!,
  set_objective_bounds_and_magnitude!,
  set_branching_rules!,
  set_branching_exp!,
  set_corecuts_cb!,
  set_facultativecuts_cb!,
  getdisaggregatedvalueofvariable
  
import ...BlockDecompositionExtras.BlockSolverExtraInterface:
  adddynamicvariable!,
  adddynamicconstraint!,
  adddynamiccut!,
  set_var_generic!,
  set_cstr_generic!,
  set_cstr_preproc!,
  set_rcsp_oracles!,
  set_rcsp_cuts!,
  set_rcsp_branching!,
  send_rcsp_to_solver!,
  addtermtodynamiccstr!,
  addtermstodynamiccstr!,
  addrhstodynamiccstr!,
  addentrytodynamiccol!,
  setobjvaluetodynamiccol!,
  getcurcostofdynvar,
  getdualofdyncstr,
  getvalueofdynvar
  
mutable struct BaPCodSolver <: AbstractMathProgSolver
  options
  vars_decomposition
  cstrs_decomposition
  sp_ids
  sp_mult
  sp_prio
  obj_magnitude
  obj_lb
  obj_ub
  vars_branching_priorities
  branching_rules
  branching_exp
  oracles
  genvarcallbacks
  gencstrcallbacks
  preproc_used_cstrs
  rcsp_oracles
  rcsp_generic_cuts
  rcsp_branching
  initial_sol
  corecuts_cb_func
  facultativecuts_cb_func
end
export BaPCodSolver

function BaPCodSolver(;kwargs...)
  include(depsfile)
  BaPCodSolver(kwargs, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

mutable struct BcMathProgModel <: AbstractMathProgModel
  inner::BcModel
  solver::BaPCodSolver
  lazycb # TODO : remove
  usercb
end

BcMathProgModel(inner, solver) =
  BcMathProgModel(inner, solver, nothing, nothing)

BcMathProgModel(s) = BcMathProgModel(BcModel(;s.options...), s)

LinearQuadraticModel(s::BaPCodSolver) = BcMathProgModel(s)
ConicModel(s::BaPCodSolver) = LPQPtoConicBridge(LinearQuadraticModel(s))
supportedcones(s::BaPCodSolver) = [:Free,:Zero,:NonNeg,:NonPos]

nb_vars(model::BcMathProgModel) = nb_vars(model.inner)
nb_cstrs(model::BcMathProgModel) = nb_cstrs(model.inner)
#getblocksolution(m::BcMathProgModel) = getblocksolution(m.inner)
getdisaggregatedvalueofvariable(model::BcMathProgModel, varidx) =
  getdisaggregatedvalueofvariable(model.inner, varidx)

################################################################################
######################   BlockDecomposition.BlockSolverInterface ########################
################################################################################
function set_constrs_decomposition!(s::BaPCodSolver, cstrs_decomposition_list)
  s.cstrs_decomposition = cstrs_decomposition_list
end

function set_vars_decomposition!(s::BaPCodSolver, vars_decomposition_list)
  s.vars_decomposition = vars_decomposition_list
end

function set_sp_mult!(s::BaPCodSolver, sp_multiplities)
  s.sp_mult = sp_multiplities
end

function set_var_branching_prio!(s::BaPCodSolver, var_priorities)
  s.vars_branching_priorities = var_priorities
end

function set_sp_ids!(s::BaPCodSolver, sp_ids)
  s.sp_ids = sp_ids
end

function set_sp_prio!(s::BaPCodSolver, sp_priorities)
  s.sp_prio = sp_priorities
end

function set_var_generic!(s::BaPCodSolver, genvarsfcts)
  s.genvarcallbacks = genvarsfcts
end

function set_cstr_generic!(s::BaPCodSolver, gencstrcbs)
  s.gencstrcallbacks = gencstrcbs
end

function set_cstr_preproc!(s::BaPCodSolver, cstrpreproc)
  s.preproc_used_cstrs = cstrpreproc
end

function set_objective_bounds_and_magnitude!(s::BaPCodSolver, magnitude, lb, ub)
  s.obj_magnitude = magnitude
  s.obj_lb = (lb > -Inf) ? lb : NaN
  s.obj_ub = (ub < Inf) ? ub : NaN
end

function set_initial_sol!(s::BaPCodSolver, solution)
    s.initial_sol = solution.initialsolutions
end

function set_branching_rules!(s::BaPCodSolver, branching_rules)
  s.branching_rules = branching_rules
end

function set_branching_exp!(s::BaPCodSolver, branching_exp)
    s.branching_exp = branching_exp
end

function set_corecuts_cb!(s::BaPCodSolver, corecutscb)
  s.corecuts_cb_func = corecutscb
end

function set_facultativecuts_cb!(s::BaPCodSolver, facultativecutscb)
  s.facultativecuts_cb_func = facultativecutscb
end

set_rcsp_oracles!(s::BaPCodSolver, oracles) = s.rcsp_oracles = oracles
set_rcsp_cuts!(s::BaPCodSolver, cuts) = s.rcsp_generic_cuts = cuts
set_rcsp_branching!(s::BaPCodSolver, branching) = s.rcsp_branching = branching

getcurrentcost(m::BcMathProgModel, varidx) = getcurrentcost(m.inner, varidx)

getcurrentub(m::BcMathProgModel, varidx) = getcurrentub(m.inner, varidx)

getcurrentlb(m::BcMathProgModel, varidx) = getcurrentlb(m.inner, varidx)

getmasterdual(m::BcMathProgModel, constridx) = getmasterdual(m.inner, constridx)

getcurcostofdynvar(m::BcMathProgModel, varname, id, sptype, spid) =
        getcurrentcost(m.inner, varname, id, sptype, spid)

getdualofdyncstr(m::BcMathProgModel, cstrname, id, sptype, spid) =
        getmasterdual(m.inner, cstrname, id, sptype, spid)

## Objective value bounds & magnitude
setobjectivevaluelb!(m::BcMathProgModel, lb) = c_set_objlb(m.inner, lb)

setobjectivevalueub!(m::BcMathProgModel, ub) = c_set_objub(m.inner, ub)

setobjectivevaluemagnitude!(m::BcMathProgModel, magnitude) =
  c_set_objmagnitude(m.inner, magnitude)

getvalueofdynvar(m::BcMathProgModel, colid::Integer) = getvalueofvariable(m.inner, colid)

################################################################################
##################### MathProgBase.SolverInterface #############################
## Loads problem data from the given file
function loadproblem!(m::BcMathProgModel, filename::AbstractString)
  println("loadproblem!(error) not implemented")
end

## Loads the provided problem data to set up the linear programming problem
##  loadproblem!(m::BcMathProgModel, A, l, u, c, lb, ub, sense)
##  which is the following LP
##        minₓ cᵀx
##        s.c. lb ≤ Ax ≤ ub
##              l ≤  x ≤ u
## sense specifies the direction of optimization problem (:Min or :Max)
function loadproblem!(m::BcMathProgModel, A, l, u, c, lb, ub, sense)
  (sense == :Max) && error("BaPCod does not support maximisation. Hint: Minimize the opposite objective")
  nbrows, nbcols = size(A)
  # Init BaPCod model
  c_init_model!(m.inner, nbrows, nbcols)

  # Objective bounds & magnitude
  !isnan(m.solver.obj_magnitude) && setobjectivevaluemagnitude!(m, m.solver.obj_magnitude)
  !isnan(m.solver.obj_lb) && setobjectivevaluelb!(m, m.solver.obj_lb)
  !isnan(m.solver.obj_ub) && setobjectivevalueub!(m, m.solver.obj_ub)

  # Registering subproblems
  (m.solver.sp_ids != nothing) && c_register_subproblems(m.inner, m.solver.sp_ids)

  # Registering variables & Constraints
  if m.solver.vars_decomposition != nothing && m.solver.cstrs_decomposition != nothing
    c_register_vars(m.inner, A, l, u, c, m.solver.vars_decomposition)
    c_register_cstrs(m.inner, A, lb, ub, m.solver.cstrs_decomposition)
  else
    c_register_vars(m.inner, A, l, u, c)
    c_register_cstrs(m.inner, A, lb, ub)
  end

  # Set subproblems multplicities
  (m.solver.sp_mult != nothing) && c_set_sp_multiplicities(m.inner, m.solver.sp_mult)
  # Set subproblem priorities
  # TODO #(m.solver.sp_prio != nothing) && c_set_sp_priorities(m.inner, m.solver.sp_prio). When doing so, check set_subproblem_priority! in wrapper>vars.jl
  #called in bc_vars.jl in function set_subproblem_priority!
  # Set variables branching priorities
  (m.solver.vars_branching_priorities != nothing) && c_vars_branching_priorities(m.inner, m.solver.vars_branching_priorities)
  # Constraints used in preproccesing
  (m.solver.preproc_used_cstrs != nothing) && c_cstrs_used_in_preproc(m.inner, m.solver.preproc_used_cstrs)
  # Branching rules
  (m.solver.branching_rules != nothing) && c_branching_rules(m.inner, m.solver.branching_rules)
  (m.solver.branching_exp != nothing) && c_branching_exp(m.inner, m.solver.branching_exp)
  # Provide initial solution
  if m.solver.initial_sol != nothing
      c_set_initial_sol_DW(m.inner,m.solver.initial_sol)
  end

    # RCSP oracles & generic cuts
    (m.solver.rcsp_oracles != nothing) && register_rcsp_oracles!(m)
    (m.solver.rcsp_generic_cuts != nothing) && register_rcsp_genericcuts!(m)
    (m.solver.rcsp_branching != nothing) && register_rcsp_branching!(m)
    # Register subproblem oracles
    (m.solver.oracles != nothing) && register_oracles!(m)
end

## Writes the current problem data to the given file.
function writeproblem(m::BcMathProgModel, filename::AbstractString)
  error("writeproblem not implemented")
end

## Returns a vector containing the lower bounds l on the variables.
getvarLB(m::BcMathProgModel) = c_get_var_lb(m.inner)

## Sets the lower bounds on the variables.
setvarLB!(m::BcMathProgModel, l) = c_set_var_lb!(m.inner, l)

##  Returns a vector containing the upper bounds u on the variables.
getvarUB(m::BcMathProgModel) = c_get_var_ub(m.inner);

## Sets the upper bounds on the variables.
setvarUB!(m::BcMathProgModel, u) = c_set_var_ub!(m.inner, u)

## Returns a vector containing the lower bounds lb on the linear constraints.
function getconstrLB(m::BcMathProgModel)
  error("getconstrLB not implemented")
end

## Sets the lower bounds on the linear constraints.
function setconstrLB!(m::BcMathProgModel, lb)
  error("setconstrLB! not implemented")
end

## Returns a vector containing the upper bounds ub on the linear constraints.
function getconstrUB(m::BcMathProgModel)
  error("getconstrUB not implemented")
end

## Sets the upper bounds on the linear constraints.
function setconstrUB!(m::BcMathProgModel, ub)
  error("setconstrUB! not implemented")
end

## Returns a vector containing the linear objective coefficients c.
getobj(m::BcMathProgModel) = c_get_var_obj(m.inner)

## Sets the linear objective coefficients.
setobj!(m::BcMathProgModel, c) = c_set_var_costs!(m.inner, c)

## Returns the full linear constraint matrix A, typically as a SparseMatrixCSC.
function getconstrmatrix(m::BcMathProgModel)
  error("getconstrmatrix! not implemented")
end

## Adds a new variable to the model, with lower bound l (-Inf if none), upper
## bound u (Inf if none), and objective coefficient objcoef. Constraint
## coefficients for this new variable are specified in a sparse format: the
## constrcoef vector contains the nonzero coefficients, and the constridx vector
## contains the indices of the corresponding linear constraints.
function addvar!(m::BcMathProgModel, constridx, constrcoef, l, u, objcoef)
  error("addvar! not implemented")
end

## Adds a new variable to the model, with lower bound l (-Inf if none),
## upper bound u (Inf if none), and objective coefficient objcoef. This is
## equivalent to calling the above method with empty arrays for the constraint
## coefficients.
function addvar!(m::BcMathProgModel, l, u, objcoef)
  error("addvar! not implemented")
end

## Adds a new linear constraint to the model, with lower bound lb (-Inf if none)
## and upper bound ub (Inf if none). Coefficients for this new constraint are
## specified in a sparse format: the coef vector contains the nonzero
## coefficients, and the varidx vector contains the indices of the corresponding
## variables.
# function addconstr!(m::BcMathProgModel, varidx, coef, lb, ub)
#   println("addconstr! not implemented")
# end
function addconstr!(m::BcMathProgModel, varidx, coef, lb, ub)
  error("addconstr! not implemented")
end

## Returns the number of linear constraints in the model.
function numlinconstr(m::BcMathProgModel)
  error("numlinconstr not implemented")
end

## Returns a vector containing the values of the linear constraints at the
## solution. This is the vector Ax.
function getconstrsolution(m::BcMathProgModel)
  error("getconstrsolution not implemented")
end

## Returns the dual solution vector corresponding to the variable bounds,
## known as the reduced costs. Not available when integer variables are present.
function getreducedcosts(m::BcMathProgModel)
  error("getreducedcosts not implemented")
end

## Returns the dual solution vector corresponding to the linear constraints.
## Not available when integer variables are present.
function getconstrduals(m::BcMathProgModel)
  error("getconstrduals not implemented")
end

## Returns the cumulative number of simplex iterations during the optimization
## process. In particular, for a MIP it returns the total simplex iterations for
## all nodes.
getsimplexiter(m::BcMathProgModel) = getstatisticcounter(m.inner, :bcCountCg)

# Provide an initial solution v to the solver, as supported. To leave values
# undefined, set them to NaN. MIP solvers should ignore provided solutions that
# are infeasible or cannot be completed to a feasible solution. Nonlinear
# solvers may use provided solutions as starting points even if infeasible.
function setwarmstart!(m::BcMathProgModel, v)
  error("setwarmstart not implemented")
end

getinfeasibilityray(m::BcMathProgModel) =
  error("Not supported. BaPCod does not return the \"Farkas\" proof of infeasibility.")

getbasis(m::BcMathProgModel) =
  error("Not supported. BaPCod does not return the basis set for the optimal solution.")

getunboundedray(m::BcMathProgModel) =
  error("Not supported. BaPCod does not return an unbounded ray of the problem.")

getbarrieriter(m::BcMathProgModel) =
  error("Not supported. BaPCod does not return the number of barrier iterations")

######################## SOLVE INTERFACE ######################################

## Returns the solution vector found by the solver.
getsolution(m::BcMathProgModel) = getaggregatedsolution(m.inner)

## Returns the objective value of the solution found by the solver. In
## particular, this may be the objective value for the best feasible solution if
## optimality is not attained or proven.
getobjval(m::BcMathProgModel) = getsolutioncost(m.inner.solution)

## Solves the optimization problem.
function optimize!(m::BcMathProgModel)
    # Register generic variable functors & constraints functors
    (m.solver.genvarcallbacks != nothing) && register_genvarcb!(m)
    (m.solver.gencstrcallbacks != nothing) && register_gencstrcb!(m)
    # Register core cuts callbacks
    (m.solver.facultativecuts_cb_func != nothing) && register_facultative_cb!(m)
    # Register facultative cuts callbacks
    (m.solver.corecuts_cb_func != nothing) && register_core_cb!(m)
    # Register lazy callbacks
    (m.lazycb != nothing) && register_lazycb!(m) # TODO : remove
    # Register user callback
    (m.usercb != nothing) && register_usercb!(m) # TODO : remove
  # Solve the model
  c_optimize(m.inner)
end

## Returns the termination status after solving. Possible values include
## :Optimal, :Infeasible, :Unbounded, :UserLimit (iteration limit or timeout),
## and :error. Solvers may return other statuses, for example, when presolve
## indicates that the model is either infeasible or unbounded, but did not
## determine which.
function status(m::BcMathProgModel)
  bestDb = getstatisticvalue(m.inner, :bcRecBestDb)
  bestInc = getstatisticvalue(m.inner, :bcRecBestInc)
  ϵ = getoptimalitygaptolerance(m.inner)
  if bestDb - ϵ <= bestInc <= bestDb + ϵ
    return :Optimal
  end
  return :UserLimit
end

## Returns the best known bound on the optimal objective value. This is used,
## for example, when a branch-and-bound method is stopped before finishing.
getobjbound(m::AbstractMathProgModel) = getstatisticvalue(m.inner, :bcRecBestDb)

## Returns the final relative optimality gap as optimization terminated.
function getobjgap(m::AbstractMathProgModel)
  error("getobjgap not implemented")
end

## Returns an object that may be used to access a solver-specific API for this
## model.
function getrawsolver(m::AbstractMathProgModel)
  error("getrawsolver not implemented")
end

##  Returns the total elapsed solution time as reported by the solver.
getsolvetime(m::AbstractMathProgModel) = getstatistictimer(m.inner, :bcTimeBaP)

##  Returns the number of nodes as reported by the solver.
getnodecount(m::AbstractMathProgModel) = getstatisticcounter(m.inner, :bcCountNodeProc)

## Sets the optimization sense of the model. Accepted values are :Min and :Max.
setsense!(m::BcMathProgModel, sense) = c_set_obj_sense!(m.inner, sense)

## Returns the optimization sense of the model.
getsense(m::BcMathProgModel) = c_get_obj_sense(m.inner)

## Returns the number of variables in the model.
function numvar(m::BcMathProgModel)
  error("numvar not implemented")
end

## Returns the total number of constraints in the model.
function numconstr(m::BcMathProgModel)
  error("numconstr not implemented")
end

## Release any resources and memory used by the model. Note that the Julia
## garbage collector takes care of this automatically, but automatic collection
## cannot always be forced. This method is useful for more precise control of
## resources, especially in the case of commercial solvers with licensing
## restrictions on the number of concurrent runs. Users must discard the model
## object after this method is invoked.
function freemodel!(m::AbstractMathProgModel)
  error("freemodel! not implemented")
end

## Sets the types of the variables to those indicated by the vector v. Valid
## types are :Int for integer, :Cont for continuous, :Bin for binary, :SemiCont
## for semicontinuous, and :SemiInt for semi-integer.
setvartype!(m::BcMathProgModel, typ::Vector{Symbol}) = c_set_var_type!(m.inner, typ)

## Returns a vector indicating the types of each variable, with values described
## above.
function getvartype(m::BcMathProgModel)
  error("getvartype not implemented")
end

## It is the philosophy of MathProgBase to not abstract over most solver
## parameters, because very few are universal, and we do not want to make it
## difficult for users who already familiar with the parameters of their solver
## of choice. However, in certain situations an abstraction over parameters is
## needed, for example, in meta-solvers which must adjust parameters of a
## subsolver inside a loop. Parameters set using these methods should override
## any parameters provided in a solver-dependent way. Solvers/models may chose
## to implement the following method:

## Sets solver-independent parameters via keyword arguments. Curent valid
## parameters are TimeLimit and Silent. TimeLimit (Float64): If the solve is not
## completed to optimality tolerances within the given number of seconds, the
## solver should return immediately with status :UserLimit. Silent (Bool): If
## set to true, the solver should disable any output to the console. If set to
## false, the parameter has no effect.
function setparameters!(m::Union{AbstractMathProgSolver, AbstractMathProgModel}; kwargs...)
    error("setparameters! not yet implemented")
end

## If these parameter-setting methods are called on an AbstractMathProgSolver,
## then they should apply to all new models created from the solver
## (but not existing models). If they are called on an AbstractMathProgModel,
## they should apply to that model only. Unrecognized parameters
## (those not listed above) should be ignored or trigger a warning message.
