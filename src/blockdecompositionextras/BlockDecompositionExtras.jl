module BlockDecompositionExtras

using ..BlockDecomposition
using JuMP

import MathProgBase
import ..BlockDecomposition.send_extras_to_solver!

import JuMP.getvalue

import ..BlockDecomposition.name,
       ..BlockDecomposition.jumpmodel

include("BlockSolverExtraInterface.jl")
using .BlockSolverExtraInterface

export  DynVariable,
        DynConstraint,
        DynVariableDict,
        DynConstraintDict,
        GenCallback,
        addvariable,
        addconstraint,
        setvarasgeneric,
        setcstrasgeneric,
        @dynconstraint,
        @dynvariable,
        @augmentconstraint,
        @augmentvariable,
        addtermtoconstraint!,
        addtermstoconstraint!,
        addrhstoconstraint!,
        addentrytocol!,
        setobjvalue!,
        isgenerated,
        getdyncurcost,
        getdyncurdual,
        constraintusedinpreprocessing!,
        getcurcostofdynvar,
        getdualofdyncstr,
        getvalue,
        name,
        jumpmodel



# export send_extras_to_solver!
function send_extras_to_solver!(model::JuMP.Model)
  send_to_solver!(model, set_var_generic!, :generic_vars)
  send_to_solver!(model, set_cstr_generic!, :generic_cstrs)
  send_to_solver!(model, set_cstr_preproc!, :cstrs_preproc)
  send_to_solver2!(model)
end

function send_to_solver2!(model)
    if applicable(send_rcsp_to_solver!, model)
       send_rcsp_to_solver!(model)
    end
end

function send_to_solver!(model::JuMP.Model, f::Function, k::Symbol)
  if haskey(model.ext, k)
    if applicable(f, model.solver, model.ext[k])
      f(model.solver, model.ext[k])
    end
  else
    error("Key $(k) does not exist in model.ext.")
  end
end

struct DynVariable
  m::JuMP.Model
  col::Int
end
struct DynConstraint
  m::JuMP.Model
  idx::Int
end
mutable struct DynVariableDict <: JuMP.JuMPContainer{JuMP.Variable, 0}
    jvar::JuMP.Variable
    vardict::Dict{Tuple, DynVariable} # Array of generated index
end
mutable struct DynConstraintDict
  name::String
  jcstr::JuMP.ConstraintRef
  cstrdict::Dict{Tuple, DynConstraint}
end

include("bjeutils.jl")
include("bjecstrs.jl")
include("bjedyncstrsvars.jl")

end # module
