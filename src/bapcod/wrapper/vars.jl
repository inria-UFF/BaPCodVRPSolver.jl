function get_var_costs!(modelptr::Ptr{Cvoid}, c::Array, nbvars::Integer)
    @bcm_ccall("getVarCosts", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                          modelptr, c, nbvars)
end

function get_var_lb!(modelptr::Ptr{Cvoid}, lb::Array, nbvars::Integer)
    @bcm_ccall("getVarLb", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                       modelptr, lb, nbvars)
end

function get_var_ub!(modelptr::Ptr{Cvoid}, ub::Array, nbvars::Integer)
    @bcm_ccall("getVarUb", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                       modelptr, ub, nbvars)
end

function get_var_type!(modelptr::Ptr{Cvoid}, t::Array, nbvars::Integer)
    @bcm_ccall("getVarType", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
                                         modelptr, t, nbvars)
end

function init_vars!(mptr::Ptr{Cvoid}, l::Array, u::Array, c::Array)
    @bcm_ccall("initVars", Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                              mptr, l, u, c)
end

function register_generic_var!(mptr::Ptr{Cvoid}, name::Symbol, column_id::Integer)
    @bcm_ccall("registerGenericVar", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
                                        mptr, name, column_id)
end

function register_var!(mptr::Ptr{Cvoid}, name::Symbol, column_id::Integer, sp_type::Integer, sp_mid::Array, var_mid::Array)
    @bcm_ccall("registerVar", Cint, ( Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint} ),
                                        mptr, name, column_id, sp_type, sp_mid, var_mid)
end

function set_subproblem_priority!(mptr::Ptr{Cvoid}, sp_bcid::Vector{Cint}, priority::Cdouble)
    @bcm_ccall("setSubproblemPriority", Cint, (Ptr{Cvoid}, Ptr{Cint}, Cdouble),
                                                    mptr, sp_bcid, priority)
end

function set_var_costs!(modelptr::Ptr{Cvoid}, c::Array, size::Integer)
    @bcm_ccall("setVarCosts", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                          modelptr, c, size)
end

# this function was not tested, neither used
function set_var_lb!(modelptr::Ptr{Cvoid}, l::Array, size::Integer)
    @bcm_ccall("setVarLb", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                       modelptr, l, size)
end

function set_var_priority_in_master!(modelptr::Ptr{Cvoid}, varname::Symbol, sp_bctype::Integer, sp_bcid::Array, priority::Cdouble)
    @bcm_ccall("setVarPriorityInMaster", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Cdouble),
                                             modelptr, varname, sp_bctype, sp_bcid, priority)
end

function set_var_priority_in_sp!(modelptr::Ptr{Cvoid}, varname::Symbol, sp_bctype::Integer, sp_bcid::Array, priority::Cdouble)
    @bcm_ccall("setVarPriorityInSp", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Cdouble),
                                            modelptr, varname, sp_bctype, sp_bcid, priority)
end

function set_var_type!(modelptr::Ptr{Cvoid}, bctyp::String, size::Integer)
    @bcm_ccall("setVarType", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
                                         modelptr, bctyp, size)
end

function set_var_ub!(modelptr::Ptr{Cvoid}, u::Array, size::Integer)
    @bcm_ccall("setVarUb", Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                       modelptr, u, size)
end

function sub_problem_mult!(mptr::Ptr{Cvoid}, mult_lb::Cint, mult_ub::Cint, sp_type::Integer, sp_bcid::Vector{Cint})
    @bcm_ccall("subProblemMult", Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cint}),
                                             mptr, mult_lb, mult_ub, sp_type, sp_bcid)
end
