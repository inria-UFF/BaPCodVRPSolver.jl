if !haskey(ENV, "BAPCOD_RCSP_LIB")
  error("The env. variable BAPCOD_RCSP_LIB was not set with the BaPCod+RCSP library path!")
end

libdevjulia = ENV["BAPCOD_RCSP_LIB"]

## ccall redef
macro bcm_ccall(func, args...)
  f = "bcInterfaceModel_$(func)"
  args = map(esc, args)
  return quote
    ccall(($f, $libdevjulia), $(args...))
  end
end

macro bcs_ccall(func, args...)
  f = "bcInterfaceSolve_$(func)"
  args = map(esc, args)
  return quote
    ccall(($f, $libdevjulia), $(args...))
  end
end

macro bcsol_ccall(func, args...)
  f = "bcSolution_$(func)"
  args = map(esc, args)
  return quote
    ccall(($f, $libdevjulia), $(args...))
  end
end

macro bcr_ccall(func, args...)
  f = "bcRCSP_$(func)"
  args = map(esc, args)
  return quote
    ccall(($f, $libdevjulia), $(args...))
  end
end

macro pisinger_knp_ccall(func, args...)
  # support removed
  # args = map(esc, args)
  # quote
  #   ccall(($func, $lib_pisinger_knp), $(args...))
  # end
end

isa_jumparray(contnr) = isa(contnr, JuMP.JuMPArray)
isa_jumpdict(contnr) = isa(contnr, JuMP.JuMPDict)
isa_jumpcontnr(contnr) = (isa_jumpdict(contnr) || isa_jumparray(contnr))
isa_array(contnr) = isa(contnr, Array)
isa_jumpcstr(contnr) = isa(contnr, JuMP.ConstraintRef)
isa_jumpvar(contnr) = isa(contnr, JuMP.Variable)

function increment_key_pos(array, key, dim)
  (key == size(array, dim)) && return increment_key_pos(array, key, dim+1)
  increment_key[dim] += 1
end


## Index to Multi Index conversion
function toArray(a)
  if isa(a, Integer)
    return [a]
  end
  if applicable(iterate, a) # we check if a is iterable
    arr = Vector{Int}()
    for i in a
      arr = vcat(arr, toArray(i))
    end
    return arr
  end
  error("Object not supported: $a of type $(typeof(a)). You can only use Integer or Iterable of supported objects (Namely, Iterable of Iterable of ... of Integer) .")
end

function sptype_to_int(sp_type::Symbol)
  if sp_type == :MIP
    return 0
  elseif sp_type == :DW_MASTER || sp_type == :B_MASTER
    return 1
  elseif sp_type == :DW_SP
    return 2
  elseif sp_type == :B_SP
    return 3
  elseif sp_type == :ALL
    return 4
  else
    error("Cannot recognize problem type : $sp_type. It must be :DW_MASTER or :DW_SP for Dantzig-Wolfe decomposition and :B_MASTER or :B_SP for Benders decomposition.")
  end
end


function createMultiIndex(array_mid::Vector{Cint}, array::Vector{Int})
  length_array = length(array)
  length_array_mid = length(array_mid)
  (length(array) > 8) && error("BaPCod does not support multi-index with more than 8 indices.")
  for i in 1:length_array_mid
    if i <= length_array
      array_mid[i] = array[i]
    else
      array_mid[i] = (i == 1) ? 0 : -1
    end
  end
end

from_index_to_BaPCodindex(id, array_mid::Vector{Cint}) = createMultiIndex(array_mid, toArray(id))
function from_index_to_BaPCodindex(id)
  bcid = Vector{Cint}(undef, 8)
  from_index_to_BaPCodindex(id, bcid)
  return bcid
end

function from_BaPCodindex_to_index(bcid)
  id = Vector{Int}()
  i = 1
  while i <= 8 && bcid[i] != -1
    push!(id, bcid[i])
    i += 1
  end
  return tuple(id...)
end
# Statistics available
# Timers
bcvalues = [:bcAverageDualSolSize,
:bcAveragePrimalSolSize,
:bcAverageSepSpBendersCutDensity,
:bcAverageSepSpDualSolDensity,
:bcFailToSolveModel,
:bcMaxDualSpaceSize,
:bcMaxIterCg,
:bcMaxNbColInMastLp,
:bcMaxNbRowInMastLp,
:bcMaxPUpd,
:bcMaxPrimalSpaceSize,
:bcMaxSepSpBendersCutDensity,
:bcMaxSepSpDualSolDensity,
:bcRecBestDb,
:bcRecBestInc,
:bcRecRootDb,
:bcRecRootInc,
:bcRecRootLpVal,
:bcSpecDbUpd,
:bcSpecIncUpd]
export bcvalues
# Values
bctimers = [
:bcTimeBaP,
:bcTimeCgSpOracle,
:bcTimeColGen,
:bcTimeMain,
:bcTimeMastMPsol,
:bcTimeMastPrepareProbConfig,
:bcTimeRedCostFixAndEnum,
:bcTimeRoot,
:bcTimeRootEval,
:bcTimeSetMast,
:bcTimeSpMPsol,
:bcTimeSpSol,
:bcTimeSpUpdateProb,
:bcTimeCutSeparation,
:bcTimeAddCutToMaster,
:bcTimeEnumMPsol,
:bcTimeSBphase1,
:bcTimeSBphase2,
:bcTimePrimalHeur]
export bctimers
# Counters
bccounters = [
:bcCountCg,
:bcCountCgSpSolverCall,
:bcCountCgT1,
:bcCountCol,
:bcCountMastSol,
:bcCountNodeGen,
:bcCountNodeProc,
:bcCountNodeTreat,
:bcCountPrimalSolSize,
:bcCountSpSol,
:bcCountCutInMaster]
export bccounters
