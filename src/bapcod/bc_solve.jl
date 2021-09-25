function c_optimize(model::BcModel)
    c_optimize(model.ptr, model.solution)
  end
  
  function c_set_initial_sol_DW(model::BcModel,solutions)
    for sol in solutions
        sol[1] = map(x -> Cint(x.col-1), sol[1])
        sol[2] = map(x -> Cdouble(x), sol[2])
  
        c_set_initial_sol_DW(model.ptr, convert(Array{Cint}, sol[1]),
                              convert(Array{Cdouble},sol[2]),length(sol[1]), sol[3][1])
    end
  
    c_provideInitialSol(model.ptr)
  end
  
  # Get the current upper bound of a variable
  function getcurrentub(model::BcModel, varidx::Integer)
    cc = Ref{Cdouble}(0.0)
    status = c_getVarCurUB(model.ptr, Cint(varidx-1), cc)
    (status != 1) && error("Cannot retrieve the current cost of the variable with varidx = $varidx.")
    return cc[]
  end
  
  # Get the current lower bound of a variable
  function getcurrentlb(model::BcModel, varidx::Integer)
    cc = Ref{Cdouble}(0.0)
    status = c_getVarCurLB(model.ptr, Cint(varidx-1), cc)
    (status != 1) && error("Cannot retrieve the current cost of the variable with varidx = $varidx.")
    return cc[]
  end
  
  # Get the dual value of constraint in the master problem
  function getmasterdual(model::BcModel, constridx::Integer)
    dd = Ref{Cdouble}(0.0)
    status = c_getMasterDual(model.ptr, Cint(constridx-1), dd)
    (status != 1) && error("Cannot retrieve the master dual of the constraint with constridx = $constridx.")
    return dd[]
  end
  
  function getmasterdual(model::BcModel, cstrname::String, cstrid, sptype::Symbol, spid)
    dd = Ref{Cdouble}(0.0)
    cstrbcid = Vector{Cint}(8)
    from_index_to_BaPCodindex(cstrid, cstrbcid)
    spbcid = Vector{Cint}(8)
    from_index_to_BaPCodindex(spid, spbcid)
    spbctype = Cint(sptype_to_int(sptype))
    status = c_getDynCstrDual(model.ptr, cstrname, cstrbcid,spbctype, spbcid, dd)
    return dd[]
  end
  
  # Get the current cost of a variable
  function getcurrentcost(model::BcModel, varidx::Integer)
    cc = Ref{Cdouble}(0.0)
    status = c_getVarCurCost(model.ptr, Cint(varidx-1), cc)
    (status != 1) && error("Cannot retrieve the current cost of the variable with varidx = $varidx.")
    return cc[]
  end
  
  function getcurrentcost(model::BcModel, varname::String, varid, sptype::Symbol, spid)
    cc = Ref{Cdouble}(0.0)
    varbcid = Vector{Cint}(8)
    from_index_to_BaPCodindex(varid, varbcid)
    spbcid = Vector{Cint}(8)
    from_index_to_BaPCodindex(spid, spbcid)
    spbctype = Cint(sptype_to_int(sptype))
    status = c_getDynVarCurCost(model.ptr, varname, varbcid,spbctype, spbcid, cc)
    return cc[]
  end
  
  # Retrieve the aggregated solution
  function getaggregatedsolution(model::BcModel)
    nbvars = nb_vars(model)
    solution = Array{Cdouble}(undef, nbvars)
    getaggregatedsolution(model, model.solution, solution)
    solution
  end
  
  function getaggregatedsolution(model::BcModel, bcsol::Ptr{Cvoid}, solarr::Array{Cdouble})
    status = c_getAggSolution(model.ptr, bcsol,solarr, length(solarr))
  #(status != 1) && @warn("Cannot retrieve the solution.")
  end
  
  # Build the solution stored in a matrix. Each row is a column.
  function getblocksolution(model::BcModel)
    nbvars = nb_vars(model)
    # For each variable we create an array containing the value of the variable
    # for each solution
    datasol = Array{Array}(undef, 0)
    bapcodsol = model.solution
    status = 1
    while status == 1
      vector = Array{Cdouble}(undef, nbvars)
      c_getValues(model.ptr, model.solution, vector, nbvars)
      for i in 1:getsolutionmultiplicity(bapcodsol)
        push!(datasol, vector)
      end
      status = getnextsolution(bapcodsol)
    end
    # Build the matrix solution with data retrieved
    solution = hcat([column for column in datasol]...)
    solution
  end
  
  function getvalueofvariable(model::BcModel, varidx::Integer)
    varsolution = 0
    bapcodsol = model.solution
    getstartsolution(model, bapcodsol)
    status = 1
    while status == 1
      value = Ref{Cdouble}(0.0)
      c_getValueOfVar(model.ptr, bapcodsol, Cint(varidx-1), value)
      varsolution += value[] * getsolutionmultiplicity(bapcodsol)
      status = getnextsolution(bapcodsol)
    end
    varsolution
  end
  
  # Build the disaggragated solution of a given variable
  function getdisaggregatedvalueofvariable(model::BcModel, varidx::Integer)
    varsolution = Array{Float64}(undef, 0)
    bapcodsol = model.solution
    getstartsolution(model, bapcodsol)
    status = 1
    while status == 1
      value = Ref{Cdouble}(0.0)
      c_getValueOfVar(model.ptr, bapcodsol, Cint(varidx-1), value)
      for i in 1:getsolutionmultiplicity(bapcodsol)
        push!(varsolution, value[])
      end
      status = getnextsolution(bapcodsol)
    end
    varsolution
  end
  
  getdisaggregatednames(model::JuMP.Model) = getdisaggregatednames(model.internalModel.inner)
  function getdisaggregatednames(model::BcModel)
    names = Array{AbstractString}(undef, 0)
    bapcodsol = model.solution
    getstartsolution(model, bapcodsol)
    status = 1
    while status == 1
      name = c_getNameOfSolForm(bapcodsol)
      name = unsafe_string(name)
      for i in 1:getsolutionmultiplicity(bapcodsol)
        push!(names, "$name#$i")
      end
      status = getnextsolution(bapcodsol)
    end
    names
  end
  
  # Get the next solution in the bapcod solution list
  # The argument solution is a pointer to a BcSolution object (C++)
  function getnextsolution(solution::Ptr{Cvoid})
    c_next(solution)
  end
  
  # Get the beginning of the solution list
  # The argument solution is a pointer to a BcSolution object (C++)
  function getstartsolution(model::BcModel, solution::Ptr{Cvoid})
    c_start(solution, model.ptr)
  end
  
  # Get the multiplicity of the solution
  # The argument solution is a pointer to a BcSolution object (C++)
  function getsolutionmultiplicity(solution::Ptr{Cvoid})
    mult = Ref{Cint}(0)
    status = c_getMultiplicity(solution, mult)
  #  (status != 1) && @warn("Cannot retrieve the multiplicity of the solution")
    return mult[]
  end
  
  # Get the solution cost
  # The argument solution is a pointer to a BcSolution object (C++)
  function getsolutioncost(solution::Ptr{Cvoid})
    cost = Ref{Cdouble}(0.0)
    status = c_getCost(solution, cost)
  # (status != 1) && @warn("Cannot retrieve the cost of the solution.")
    return cost[]
  end
  
  # Get the solution true cost
  # The argument solution is a pointer to a BcSolution object (C++)
  function getsolutiontruecost(solution::Ptr{Cvoid})
    truecost = Cdouble(0)
    status = c_getTrueCost(solution, truecost)
   (status != 1) && @warn("Cannot retrieve the true cost of the solution.")
   return truecost[]
  end
  
  function getoptimalitygaptolerance(model::BcModel)
    ogt = Ref{Cdouble}(0.0)
    status = c_getOptimalityGapTolerance(model.ptr, ogt)
    (status != 1) && error("Cannot retrieve the optimality gap tolerance.")
    return ogt[]
  end
  
  getstatisticvalue(m::BcModel, key::Symbol) =
    c_getStatisticValue(m.ptr, "$key")
  
  getstatisticcounter(m::BcModel, key::Symbol) =
    c_getStatisticCounter(m.ptr, "$key")
  
  getstatistictimer(m::BcModel, key::Symbol) =
    c_getStatisticTime(m.ptr, "$key")
  
  function getstatistic(m::JuMP.Model, key::Symbol)
    bcm = internalmodel(m).inner
    (!isempty(filter(x -> x != 0, bcvalues .== key))) && return getstatisticvalue(bcm, key)
    (!isempty(filter(x -> x != 0, bccounters .== key))) && return getstatisticcounter(bcm, key)
    (!isempty(filter(x -> x != 0, bctimers .== key))) && return getstatistictimer(bcm, key)
    error("Unknown statistic $key.")
  end
  
  # Get the current upper bound of a variable
  function getenumeratedcols(model::BcModel, bcsol::Ptr{Cvoid})
    nbevals = Ref{Cint}(0)
    status = c_enumerateAllColumns(model.ptr, bcsol, nbevals)
    return nbevals[]
  end
  