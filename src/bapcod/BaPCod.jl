module BaPCod
  using JuMP
  using ..BlockDecomposition
  using ..BlockDecompositionExtras
  using SparseArrays

  import Base.convert, Base.show, Base.copy, Base.pointer
  using MathProgBase.SolverInterface
  using ..BlockDecomposition.BlockSolverInterface
  using ..BlockDecompositionExtras.BlockSolverExtraInterface

  # Heuristics
  export solvewith_pisinger_minknap,
         solvewith_pisinger_bouknap,
         set_solution_to_homogeneous_syst,
         set_solution_to_homogeneous_syst_with_eq,
         get_oracle_phaseofstageapproach,
         getsolutioncost,
         getsolutiontruecost,
         getstatistic,
         getdisaggregatednames

  export getcstrindex, getvarindex, setobjvaluetodynamiccol!

  include("bc_common.jl")

  include("wrapper/callbacks.jl")
  include("wrapper/constrs.jl")
  include("wrapper/model.jl")
  include("wrapper/oracles.jl")
  include("wrapper/rcsp.jl")
  include("wrapper/vars.jl")
  include("wrapper/solve.jl")

  include("bc_heuristics.jl")
  include("bc_model.jl")
  include("bc_vars.jl")
  include("bc_constrs.jl")
  include("BaPCodSolverInterface.jl")
  include("bc_solve.jl")
  include("bc_rcsp.jl") # done
  include("bc_oracles.jl")
  include("bc_callbacks.jl") # done
  include("bc_gencallback.jl")
  include("bc_spec_callbacks.jl") # done

end # module
