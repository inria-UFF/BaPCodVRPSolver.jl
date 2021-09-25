module BlockSolverInterface

using Base.Meta
using Compat

abstract type OracleSolverData end
export OracleSolverData

abstract type OracleCallbackData end
export OracleCallbackData

# Interfaces
"""
    set_constrs_decomposition!(s::AbstractMathProgSolver, data::Array)

sends to the solver `s` in which subproblems are the constraints.
Each element of the `data` array is a `Tuple` containing
`(constr_name::Symbol, constr_id::Tuple, sp_type::Symbol, sp_id::Tuple)`.
`constr_name` and `constr_id` are the name and the index of the constraint in the JuMP model.
`sp_type` and `sp_id` are the type and the index of the subproblem to which the constraint is assigned.
"""
set_constrs_decomposition!() = nothing
export set_constrs_decomposition!

"""
    set_vars_decomposition!(s::AbstractMathProgSolver, data::Array)

sends to the solver `s` in which subproblems are the variables.
Each element of the `data` array is a `Tuple` containing
`(var_name::Symbol, var_id::Tuple, sp_type::Symbol, sp_id::Tuple)`.
`var_name` and `var_id` are the name and the index of the variable in the JuMP model.
`sp_type` and `sp_id` are the type and the index of the subproblem to which the variable is assigned.
"""
set_vars_decomposition!() = nothing
export set_vars_decomposition!

set_sp_ids!() = nothing
export set_sp_ids!

"""
    set_sp_mult!(s::AbstractMathProgSolver, multiplicities::Array)

sends to the solver `s` the multiplicity of each subproblem. Each element of
the `multiplicities` array is a `Tuple` containing `(sp_id::Tuple, sp_type::Symbol, mult_lb, mult_ub)`
with `mult_lb` the lower bound of the multiplicity and `mult_ub` the upper
bound of the multiplicity.
"""
set_sp_mult!() = nothing
export set_sp_mult!

"""
    set_sp_prio!(s::AbstractMathProgSolver, priorities::Array)

sends to the solver `s` the list of subproblem priorities. Each element of
the `data` array is a `Tuple` containing `(sp_id::Tuple, sp_type::Symbol, sp_prio)`.
"""
set_sp_prio!() = nothing
export set_sp_prio!

"""
    set_var_branching_prio!(s::AbstractMathProgSolver, priorities::Array)

sends to the solver `s` the list of variables branching priorities.
The number stored at the row `i` is the branching priority of the variable
stored at the column `i` in the JuMP model.
"""
set_var_branching_prio!() = nothing
export set_var_branching_prio!

"""
    set_oracles!(s::AbstractMathProgSolver, oracles::Array)

sends to the solver `s` the `list` with subproblems and oracles functions.
Each element of the `oracles` array is a `Tuple` containing
`(sp_id::Tuple, sp_type::Symbol, oracle::Function)` with `oracle`, the
function defined by the user.
"""
set_oracles!() = nothing
export set_oracles!

"""
    set_objective_bounds_and_magnitude!(s::AbstractMathProgSolver, magn, lb, ub)


sends to the solver `s` the magnitude `magn`, the lower bound `lb`
and the upper bound `ub` of the objective function.
"""
set_objective_bounds_and_magnitude!() = nothing
export set_objective_bounds_and_magnitude!


"""
    set_oraclesolution_newsolution(o::OracleSolverData)

creates a new solution in the oracle solver solution. It is usefull, if the
user wants to return several solutions.
"""
set_oraclesolution_newsolution() = nothing
export set_oraclesolution_newsolution

"""
    set_oraclesolution_solution(o::OracleSolverData, x::JuMP.Variable, v::Real)

Set the value of the variable `x` to `v` in the oracle solver solution stored
in the `OracleSolverData` object `o`.
"""
set_oraclesolution_solution() = nothing
export set_oraclesolution_solution

"""
    set_oraclesolution_objval(o::OracleSolverData, v::Real)


Sets the objective value stored in the `OracleSolverData` object `o` to `v`.
"""
set_oraclesolution_objval() = nothing
export set_oraclesolution_objval

"""
    get_oracle_phaseofstageapproach(o::OracleSolverData)

Returns the phase of stage approach.
"""
get_oracle_phaseofstageapproach() = nothing
export get_oracle_phaseofstageapproach

"""
    getcurrentcost(m::AbstractMathProgModel, vcol::Integer)


returns the current cost of the `vcol` ``{}^{th}`` variable.
"""
getcurrentcost() = nothing
export getcurrentcost

"""
    getcurrentub(m::AbstractMathProgModel, vcol::Integer)


returns the current ub of the `vcol` ``{}^{th}`` variable.
"""
getcurrentub() = nothing
export getcurrentub

"""
    getcurrentlb(m::AbstractMathProgModel, vcol::Integer)


returns the current lb of the `vcol` ``{}^{th}`` variable.
"""
getcurrentlb() = nothing
export getcurrentlb

"""
    getmasterdual(m::AbstractMathProgModel, crow::Integer)


returns the current dual of the `crow` ``{}^{th}`` variable.
"""
getmasterdual() = nothing
export getmasterdual

"""
    getdisaggregatedvalueofvariable(m::AbstractMathProgModel, vcol::Integer)

returns the disaggregated value of the `vcol` ``{}^{th}`` variable.
"""
getdisaggregatedvalueofvariable() = nothing
export getdisaggregatedvalueofvariable

"""
    set_branching_rules!(s::AbstractMathProgSolver, rules::Dict{Symbol, Any})

sends to the solver `s` a Dictionary `rules`. The key is the name of the branching rule
and the content is an array of branching instances. The array contains Tuple of variables name and parameters.
Parameters are stored in an array of tuple (name_of_parameter, value_of_parameter).

For instance,

    rules = (:branching_rule_name => [(:x, [(:priority, 1)]), (:y, [(:priority, 1)])])

Names of branching rules depend on solvers.
"""
set_branching_rules!() = nothing
export set_branching_rules!

set_branching_exp!() = nothing
export set_branching_exp!

set_corecuts_cb!() = nothing
export set_corecuts_cb!

set_facultativecuts_cb!() = nothing
export set_facultativecuts_cb!

getvalvariable() = nothing
export getvalvariable


end
