module BlockDecomposition

include("BlockSolverInterface.jl")
using .BlockSolverInterface

#importall JuMP

# This issue on JuliaLang ends with this comment (Closing as "don't use Requires.jl")
# https://github.com/JuliaLang/julia/issues/20909

# using Requires
using JuMP
using MathProgBase
using SparseArrays

# DEPENDENCY ON BlockDecompositionExtras IS NOW IMPLEMENTED BY REDEFINING THE FUNCTION
# send_extras_to_solver

# @require BlockDecompositionExtras begin
#        using BlockDecompositionExtras
# end

# DEPENDENCY ON CPLEX IS MENDATORY (until we add a similar hook to be redefined in CPLEX)

# @require CPLEX begin
#       using CPLEX
# end

export  BlockModel,
        BlockIdentificationData,
        BlockDecompositionData,
        defaultBlockIdeData,
        getcurcost,
        getcurub,
        getcurlb,
        getcurdual,
        getdisaggregatedvalue,
        objectivevaluemagnitude,
        objectivevalueupperbound,
        objectivevaluelowerbound,
        branchingpriorityinmaster,
        branchingpriorityinsubproblem,
        add_Dantzig_Wolfe_decomposition,
        add_Dantzig_Wolfe_decomposition_on_variables,
        add_Benders_decomposition,
        addspmultiplicity,
        addsppriority,
        addbranching,
        show

# Oracles
export OracleSolverData,
       addoracletosp!,
       add_oracle_to_DWsp!,
       add_oracle_to_Bsp!,
       setsolutionvalue,
       setsolutionbestobjval,
       addsolution,
       getspid,
       getsptype,
       addsolution,
       addfacultativecutcallback,
       addcorecutcallback,
       getphaseofstageapproach,
       getvalue

export name, jumpmodel

import Base.convert, Base.show, Base.copy, Base.pointer

include("bjutils.jl")
include("bjreport.jl")
include("bjprint.jl")
include("bjmodel.jl")
include("bjexpand.jl")
include("bjoracles.jl")
include("bjsolve.jl")
include("bjdecomposition.jl")
#include("bjcplex.jl")
include("bjcallbacks.jl")

end # module
