function c_addToOracleSol(mptr::Ptr{Cvoid},solptr,var_id,var_val)
    @bcs_ccall("addToOracleSol", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cdouble),
                                        mptr,solptr,var_id,var_val)
end

function c_getOraclePhaseOfStage(solptr::Ptr{Cvoid})
  @bcs_ccall("getOraclePhaseOfStage", Cint, (Ptr{Cvoid},), solptr)
end

function c_initOracle!(mptr::Ptr{Cvoid}, oraclesolver_fct, sp_type, sp_bcid, m, oracle_id)
    @bcs_ccall("initOracle", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cint}, Any, Cint),
                                mptr, oraclesolver_fct, sp_type, sp_bcid, m, oracle_id)
end

function c_newOracleSol!(solptr::Ptr{Cvoid})
  @bcs_ccall("newOracleSol", Cint, (Ptr{Cvoid},), solptr)
end

function c_setObjValOfOracleSol(solptr::Ptr{Cvoid},objval)
    @bcs_ccall("setObjValOfOracleSol", Cint, (Ptr{Cvoid}, Cdouble),
                                          solptr, objval)
end

function c_setSolToHomogeneousSyst(solptr::Ptr{Cvoid})
  @bcs_ccall("setSolToHomogeneousSyst", Cint, (Ptr{Cvoid},), solptr)
end

function c_setSolToHomogeneousSystWithEq(solptr::Ptr{Cvoid})
  @bcs_ccall("setSolToHomogeneousSystWithEq", Cint, (Ptr{Cvoid},), solptr)
end
