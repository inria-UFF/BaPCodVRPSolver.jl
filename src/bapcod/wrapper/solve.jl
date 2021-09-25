function c_getAggSolution(mptr::Ptr{Cvoid}, bcsol, solarr, length_solarr)
    @bcsol_ccall("getAggSolution", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                            mptr, bcsol, solarr, length_solarr)
end

function c_getCost(solution::Ptr{Cvoid}, cost)
  @bcsol_ccall("getCost", Cint, (Ptr{Cvoid}, Ref{Cdouble}),
                                          solution, cost)
end

function c_getDynVarCurCost(mptr::Ptr{Cvoid}, varname, varbcid, spbctype, spbcid, dd)
    @bcs_ccall("getDynVarCurCost", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cint}, Cint,
                                              Ptr{Cint}, Ref{Cdouble}),
                                              mptr, varname, varbcid,
                                              spbctype, spbcid, dd)
end

function c_getDynCstrDual(mptr::Ptr{Cvoid}, cstrname, cstrbcid, spbctype, spbcid, cc)  
    @bcs_ccall("getDynCstrDual", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cint}, Cint,
                                              Ptr{Cint}, Ref{Cdouble}),
                                              mptr, cstrname, cstrbcid,
                                              spbctype, spbcid, cc)
end

function c_getMasterDual(mptr::Ptr{Cvoid}, cstrs_id, dd)
    @bcs_ccall("getMasterDual", Cint, (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                mptr, cstrs_id, dd)
end

function c_getMultiplicity(solution::Ptr{Cvoid}, mult)
  @bcsol_ccall("getMultiplicity", Cint, (Ptr{Cvoid}, Ref{Cint}),
                                        solution, mult)
end

function c_getNameOfSolForm(bapcodsol::Ptr{Cvoid})
    @bcsol_ccall("getNameOfSolForm", Cstring, (Ptr{Cvoid},),bapcodsol)
end

function c_getStatisticCounter(mptr::Ptr{Cvoid}, key::String)
  @bcs_ccall("getStatisticCounter", Clong, (Ptr{Cvoid}, Ptr{Cchar}), mptr, key)
end

function c_getStatisticTime(mptr::Ptr{Cvoid}, key::String)
  @bcs_ccall("getStatisticTime", Clong, (Ptr{Cvoid}, Ptr{Cchar}), mptr, key)
end

function c_getStatisticValue(mptr::Ptr{Cvoid}, key::String)
  @bcs_ccall("getStatisticValue", Cdouble, (Ptr{Cvoid}, Ptr{Cchar}), mptr, key)
end

function c_getTrueCost(solution::Ptr{Cvoid}, truecost::Float64)
  @bcsol_ccall("getTrueCost", Cint, (Ptr{Cvoid}, Ref{Cdouble}),
                                        solution, truecost)
end

function c_getValues(mptr::Ptr{Cvoid}, solution::Ptr{Cvoid}, vector, nbvars::Integer)
    @bcsol_ccall("getValues", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                     mptr, solution, vector, nbvars)
end

function c_getValueOfVar(mptr::Ptr{Cvoid}, bapcodsol, var_idx::Integer, value)
  @bcsol_ccall("getValueOfVar", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                        mptr, bapcodsol, var_idx, value)
end

function c_getVarCurCost(mptr::Ptr{Cvoid}, var_id, cc)
  @bcs_ccall("getVarCurCost", Cint, (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                    mptr, var_id, cc)
end

function c_getVarCurLB(mptr::Ptr{Cvoid}, var_id::Integer, cc::Ref{Cdouble})
    @bcs_ccall("getVarCurLB", Cint, (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                    mptr, var_id, cc)
end

function c_getVarCurUB(mptr::Ptr{Cvoid}, vars_idx::Integer, cc::Ref{Cdouble})
  @bcs_ccall("getVarCurUB", Cint, (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                    mptr, vars_idx, cc)
end

function c_getOptimalityGapTolerance(mptr::Ptr{Cvoid}, ogt::Ref{Cdouble})
  @bcs_ccall("getOptimalityGapTolerance", Cint, (Ptr{Cvoid}, Ref{Cdouble}), mptr, ogt)
end

function c_next(solution::Ptr{Cvoid})
    @bcsol_ccall("next", Cint, (Ptr{Cvoid},), solution)
end

function c_optimize(modelptr::Ptr{Cvoid}, solution::Ptr{Cvoid})
  @bcs_ccall("optimize", Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), modelptr, solution)
end

function c_provideInitialSol(mptr::Ptr{Cvoid})
  @bcs_ccall("provideInitialSol", Cvoid, (Ptr{Cvoid},), mptr)
end

function c_set_initial_sol_DW(mptr::Ptr{Cvoid}, vars_id, vars_val, nb_var::Integer, sp_id::Integer)
    @bcs_ccall("setInitialSol", Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble},Cint,Cint),
                                        mptr, vars_id, vars_val, nb_var, sp_id)
end

function c_start(solution::Ptr{Cvoid}, mptr)
  @bcsol_ccall("start", Cint, (Ptr{Cvoid}, Ptr{Cvoid}), solution, mptr)
end

function c_print(solution::Ptr{Cvoid})
    @bcsol_ccall("print", Cint, (Ptr{Cvoid},), solution)
end

function c_enumerateAllColumns(mptr::Ptr{Cvoid}, solution::Ptr{Cvoid}, nbcols::Ref{Cint})
  #enumerateAllColumns(InterfaceModel* m, BcSolution* s, int & nbEnumeratedCols);
  @bcsol_ccall("enumerateAllColumns", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ref{Cint}),
    mptr, solution, nbcols)
end

function c_get_sol_resource_consumption(c_sol::Ptr{Cvoid}, res_id::Integer)
    nbnodes = @bcsol_ccall("getNbNodes", Cint, (Ptr{Cvoid},), c_sol)
    nodes_id = zeros(Cint, nbnodes)
    consumption = zeros(Cdouble, nbnodes)
    status = @bcsol_ccall("getResConsumption", Cint,
                         (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint),
                         c_sol, consumption, nodes_id, Cint(nbnodes), Cint(res_id))
    (status != 1) && error("cannot retrieve resource consumption of nodes.")
    return nodes_id .+ 1, consumption
end

function c_get_sol_arcs_ids(c_sol::Ptr{Cvoid})
    nbnodes = @bcsol_ccall("getNbNodes", Cint, (Ptr{Cvoid},), c_sol)
    arcs_ids = zeros(Cint, nbnodes - 1)
    status = @bcsol_ccall("getArcsIds", Cint,
                         (Ptr{Cvoid}, Ptr{Cint}, Cint),
                         c_sol, arcs_ids, Cint(nbnodes - 1))
    (status != 1) && error("cannot retrieve arcs ids.")
    return arcs_ids 
end

function c_sol_print(c_sol::Ptr{Cvoid})
  status = @bcsol_ccall("print", Cint, (Ptr{Cvoid},), c_sol)
  return status
end

function c_get_prob_first_id(c_sol::Ptr{Cvoid})
  status = @bcsol_ccall("getProblemFirstId", Cint, (Ptr{Cvoid},), c_sol)
  return status
end

