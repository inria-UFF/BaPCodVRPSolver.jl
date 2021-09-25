function DynVariable(m::JuMP.Model)
  m.ext[:colscounter] += 1
  col = m.ext[:colscounter]
  return DynVariable(m, col)
end

function DynConstraint(m::JuMP.Model)
  m.ext[:rowscounter] += 1
  row = m.ext[:rowscounter]
  return DynConstraint(m, row)
end

DynVariableDict(v::JuMP.Variable) =
  DynVariableDict(v, Dict{Tuple, DynVariable}())

DynConstraintDict(n::String, c::JuMP.ConstraintRef) =
  DynConstraintDict(n, c, Dict{Tuple, DynConstraint}())

Base.getindex(d::DynVariableDict, t...) = d.vardict[t]
Base.setindex!(d::DynVariableDict, value::DynVariable, t...) =
  (d.vardict[t...] = value)

isgenerated(d::DynVariableDict, t...) = haskey(d.vardict, t)

Base.getindex(d::DynConstraintDict, t...) = d.cstrdict[t]
Base.setindex!(d::DynConstraintDict, value::DynConstraint, t...) =
  (d.cstrdict[t...] = value)

isgenerated(d::DynConstraintDict, t...) = haskey(d.cstrdict, t)

function setvarasgeneric(x, f::Function)
  if !isa_JuMPvar(x)
    info = """ The second argument of setvarasgeneric must be a JuMP.Variable.
               Make sure that you have declared just one variable without index.
               Example: @variable(model, name >= 0, type)
           """
    error(info, "Second argument is a $(typeof(x)).")
  end
  model = jumpmodel(x)
  if !haskey(model.ext, :generic_vars)
    info = """ The model does not support dynamic variables.
               Make sure that you use a model built with BlockDecomposition.
           """
    error(info)
  end
  model.ext[:generic_vars][Symbol(getname(x))] = (x, f)
end

function setcstrasgeneric(c, name::String, f::Function)
  if !isa_JuMPcstr(c)
    info = """ The second argument of setcstrasgeneric must be a  JuMP.ConstraintRef.
               Make sure that you have declared just one constraint without index.
               Example: @constraint(model, name, 0 == 0)
           """
    bjerror(info, "Second argument is a $(typeof(c))")
  end
  model = jumpmodel(c)
  if !haskey(model.ext, :generic_cstrs)
    info = """ The model does not support dynamic constraints.
               Make sure that you use a model built with BlockDecomposition.
           """
    error(info)
  end
  model.ext[:generic_cstrs][c.idx] = (c, name, f)
end

macro dynvariable(args...)
  nbargs = length(args)
  if length(args) < 3 || length(args) > 4
    error("Incorrect number of arguments.")
  end

  model = args[1]
  varexpr = string(args[2])
  varfunctor = args[3]
  vartype = (nbargs == 4) ? args[4] : :Continue

  # Check syntax & get the name
  varname = :genericvar
  syntax1 = match(r"^(.* <=)? (?P<varname>[a-zA-Z]*)\[\] <= .*$", varexpr)
  syntax2 = match(r"^(?P<varname>[a-zA-Z]*)\[\] >= .*$", varexpr)
  if syntax1 != nothing
    varname = Symbol(syntax1[:varname])
  elseif syntax2 != nothing
    varname = Symbol(syntax2[:varname])
  else
    error("Incorrect syntax.")
  end
  vargenerator = parse(replace(varexpr, r"\[\]", ""))
  if vartype == :Continue
    exp = quote var = @variable($model, $vargenerator) end
  else
    exp = quote var = @variable($model, $vargenerator, $vartype) end
  end
  exp = quote
    $exp
    setvarasgeneric(var, $varfunctor)
    $varname = DynVariableDict(var)
  end
  return esc(exp)
end

macro dynconstraint(args...)
  nbargs = length(args)
  if nbargs != 3
    bjerror("Incorrect number of arguments.")
  end

  model = args[1]
  cstrexpr = string(args[2])
  cstrfunctor = args[3]

  cstrname = :genericcstr
  syntax = match(r"^(?P<cstrname>[a-zA-Z]*)\[\]$", cstrexpr)
  if syntax == nothing
    error("Incorrect syntax")
  end
  cstrname = Symbol(syntax[:cstrname])
  name = String(syntax[:cstrname])
  # Check syntax & get the name
  exp = quote
    cstr = @constraint($model, $cstrname, 0 == 0)
    setcstrasgeneric(cstr, $name, $cstrfunctor)
    $cstrname = DynConstraintDict($name, cstr)
  end
  return esc(exp)
end

macro augmentvariable(args...)
  cbdata = args[1]
  vargenerator = args[2]
  dynvar = args[2].args[1]
  indexsets = args[2].args[2:end]
  err = :()
  for indexset in indexsets
    err = quote
      $err
      if !applicable(start, $indexset)
        error("\e[31m error \e[00m")
      end
    end
  end

  counter = 1
  countervars = [:a, :b, :c, :e, :d, :f, :g, :h]
  tuple = "tuple("
  for i in length(indexsets):-1:1
    tuple = "$tuple $(string(countervars[i])),"
  end
  tuple = "$tuple)"
  tuple = parse(tuple)
  exp = quote addvariable($dynvar, $tuple) end
  for indexset in reverse(indexsets)
    exp = quote
      for $(countervars[counter]) in $indexset
        $exp
      end
    end
    counter += 1
  end
  code = quote
    $err
    $exp
  end
  esc(code)
end

macro augmentconstraint(args...)
  cbdata = args[1]
  cstrgenerator = args[2]
  dyncstr = args[2].args[1]
  indexsets = args[2].args[2:end]

  counter = 1
  countervars = [:a, :b, :c, :e, :d, :f, :g, :h]
  tuple = "tuple("
  for i in length(indexsets):-1:1
    tuple = "$tuple $(string(countervars[i])),"
  end
  tuple = "$tuple)"
  tuple = Meta.parse(tuple)
  exp = quote
    addconstraint($cbdata, $dyncstr, $tuple)
  end
  for indexset in reverse(indexsets)
    exp = quote
      for $(countervars[counter]) in $indexset
        $exp
      end
    end
    counter += 1
  end
  esc(exp)
end

function addvariable(d::DynVariableDict, index::Tuple)
  model = jumpmodel(d.jvar)
  DW_dec_on_vars_f = model.ext[:block_decomposition].DantzigWolfe_decomposition_on_vars_fct
  (sp_type, sp_id) = (:DW_MASTER, (0,))
  if DW_dec_on_vars_f != nothing
    (sp_type, sp_id) = DW_dec_on_vars_f(Symbol(getname(d.jvar)), index)
  end
  if isa(sp_id, Integer) sp_id = (sp_id,) end
  dvar = DynVariable(model)
  d[index] = dvar
  colid = Cint(dvar.col)
  if !applicable(adddynamicvariable!, internalmodel(model), getname(d.jvar), index, colid, sp_type, sp_id)
     error("The solver does not support dynamic variables.")
  end
  adddynamicvariable!(internalmodel(model), getname(d.jvar), index, colid, sp_type, sp_id)
  return dvar
end

function addconstraint(cb::MathProgBase.MathProgCallbackData, d::DynConstraintDict, index::Tuple)
  model = jumpmodel(d.jcstr)
  dcstr = DynConstraint(model)
  d[index] = dcstr
  if !applicable(adddynamiccut!, internalmodel(model), d.name, index, dcstr.idx, cb.cutlist)
    error("The solver does not support dynamic constraints in cut callbacks.")
  end
  adddynamiccut!(internalmodel(model), d.name, index, dcstr.idx, cb.cutlist)
  return dcstr
end

function addconstraint(cb::BlockDecomposition.BlockSolverInterface.OracleSolverData, d::DynConstraintDict, index::Tuple)
  model = jumpmodel(d.jcstr)
  dcstr = DynConstraint(model)
  d[index] = dcstr
  if !applicable(adddynamicconstraint!, internalmodel(model), d.name, index, dcstr.idx, cb)
    error("The solver does not support dynamic constraints in the oracle.")
  end
  adddynamicconstraint!(internalmodel(model), d.name, index, dcstr.idx, cb)
  return dcstr
end

function addtermtoconstraint!(cb::MathProgBase.MathProgCallbackData, var::Union{JuMP.Variable, DynVariable}, coeff)
  colid = Cint(var.col)
  coeff = Cdouble(coeff)
  if !applicable(addtermtodynamiccstr!, cb, colid, coeff)
    error("Cannot add term to constraint")
  end
  addtermtodynamiccstr!(cb, colid, coeff)
  return
end

function addtermstoconstraint!(cb::MathProgBase.MathProgCallbackData, varcols::Vector{Int}, coeff)
  addtermstodynamiccstr!(cb, varcols, coeff) 
  return
end

function addentrytocol!(cb::MathProgBase.MathProgCallbackData, cstr::Union{JuMP.ConstraintRef, DynConstraint}, coeff)
  coeff = Cdouble(coeff)
  cstrid = Cint(cstr.idx)
  if !applicable(addentrytodynamiccol!, cb, cstrid, coeff)
    error("Cannot set entry to the dynamic column.")
  end
  addentrytodynamiccol!(cb, cstrid, coeff)
  return
end

function setobjvalue!(cb::MathProgBase.MathProgCallbackData, value)
  if !applicable(setobjvaluetodynamiccol!, cb, value)
    error("Cannot set objective value to the dynamic variable.")
  end
  setobjvaluetodynamiccol!(cb, value)
  return
end

function addrhstoconstraint!(cb::MathProgBase.MathProgCallbackData, sense, rhs)
  ssense = string(sense)
  schar = 'L'
  if ssense == "<="
    schar = 'L'
  elseif ssense == ">="
    schar = 'G'
  elseif ssense == "=="
    schar = 'E'
  else
    error("Unknown sense $ssense. Must be <=, >= or ==.")
  end
  rhs = Cdouble(rhs)
  if !applicable(addrhstodynamiccstr!, cb, schar, rhs)
    error("Cannot add rhs to constraint")
  end
  addrhstodynamiccstr!(cb, schar, rhs)
  return
end

function getdyncurcost(x::DynVariableDict, t...)
  model = jumpmodel(x.jvar)
  DW_dec_on_vars_f = model.ext[:block_decomposition].DantzigWolfe_decomposition_on_vars_fct
  if DW_dec_on_vars_f != nothing
    (sp_type, sp_id) = DW_dec_on_vars_f(Symbol(getname(x.jvar)), t)
  else
    (sp_type, sp_id) = (:DW_MASTER, (0,))
  end
  if applicable(getcurcostofdynvar, internalmodel(model), getname(x.jvar), t, sp_type, sp_id)
    getcurcostofdynvar(internalmodel(model), getname(x.jvar), t, sp_type, sp_id)
  else
    error("Cannot get the cost of the dynamic variable.")
  end
end

function getdyncurdual(x::DynConstraintDict, t...)
  model = jumpmodel(x.jcstr)
  DW_dec_f = model.ext[:block_decomposition].DantzigWolfe_decomposition_fct
  if DW_dec_f != nothing
    (sp_type, sp_id) = DW_dec_f(Symbol(x.name), t)
  else
    (sp_type, sp_id) = (:DW_MASTER, (0,))
  end
  if applicable(getdualofdyncstr, internalmodel(model), x.name, t, sp_type, sp_id)
    getdualofdyncstr(internalmodel(model), x.name, t, sp_type, sp_id)
  else
    error("Cannot get the cost of the dynamic variable.")
  end
end

function getvalue(x::DynVariable)
  model = jumpmodel(x)
  if applicable(getvalueofdynvar, internalmodel(model), x.col)
    return getvalueofdynvar(internalmodel(model), x.col)
  else
    error("Cannot retrieve the value of the dynamic variable.")
  end
end
