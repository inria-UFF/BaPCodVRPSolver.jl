function c_cstr_used_in_preprocessing!(mptr::Ptr{Cvoid}, cstr_name::Ptr{Cchar}, sp_bctype::Cint, sp_bcid::Vector{Cint}, used::Cint)
    @bcm_ccall("cstrUsedInPreprocessing", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint,
                                                      Ptr{Cint}, Cint),
                                                      m.ptr, cstr_name, sp_bctype,
                                                      sp_bcid, used)
end

function init_cstrs!(mptr::Ptr{Cvoid}, coeffs::Array, rows_id::Array, nonzeros::Array, lb::Array, ub::Array)
    @bcm_ccall("initCstrs", Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                                          mptr, coeffs, rows_id, nonzeros, lb, ub)
end

function register_cstr!(mptr::Ptr{Cvoid}, name::Symbol, row_id::Cint, sp_type::Integer, sp_mid::Array, cstr_mid::Array)
    @bcm_ccall("registerCstr", Cint, ( Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint} ),
                                        mptr, string(name), row_id, sp_type, sp_mid, cstr_mid)
end

function register_generic_cstr!(mptr::Ptr{Cvoid}, name::Symbol, row_id::Cint)
    @bcm_ccall("registerGenericCstr", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
                                                   mptr, name, row_id)
end
