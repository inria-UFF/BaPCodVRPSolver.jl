"""
    constraintusedinpreprocessing!(c::JuMP.JuMPContainer, sp::Tuple{Symbol, Union{Integer, Tuple}}, used::Bool)

Constraints in the JuMP container `c` of the subproblem `sp` will be used in the
preprocessing according to the value of `used`.

```julia
  constraintusedinpreprocessing!(cov, (:DW_MASTER, 1), false)
```
Constraints `cov` in the Danztig-Wolfe subproblem with id `1` are not used in
the preprocessing.
"""
function constraintusedinpreprocessing!(c::JuMP.JuMPContainer, sp::Tuple{Symbol, Union{Integer, Tuple}}, used::Bool)
  # get the model
  model = jumpmodel(c)
  name = lookfornameofJuMPobj(model, c)
  if !haskey(model.ext, :cstrs_preproc)
    info = """ Cannot (des)activate preproccesing on a constraint.
               Make sure that you use a model built with BlockDecomposition.
           """
    error(info)
  end
  model.ext[:cstrs_preproc][(name, sp)] = used
end
