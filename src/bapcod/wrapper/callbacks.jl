## Lazy & user callbacks
function c_get_sep_cb_type!(modelinnerptr::Ptr{Cvoid}, idcb::Cint)
    @bcs_ccall("getSepCbType", Cchar, (Ptr{Cvoid}, Cint),
                                    modelinnerptr, idcb)
end

function wbcs_init_sep_routine!(c_model::Ptr{Cvoid}, cb_c_fct::Ptr{Cvoid},
                                type_cb::Char, mpb_model::Any)
  @bcs_ccall("initSepRoutine", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Any, Cchar),
                                      c_model, cb_c_fct,  mpb_model, type_cb)
end

function wbcs_add_sep_cut!(c_model::Ptr{Cvoid}, id_cb::Cint, cutslist::Ptr{Cvoid},
                          vars_coeffs::Vector{Cdouble}, vars_idx::Vector{Cint},
                          nb_vars::Cint, sense::Char, rhs::Cdouble)
  @bcs_ccall("addSepCut", Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Ptr{Cdouble},
                                 Ptr{Cint}, Cint, Cchar, Cdouble),
                                 c_model, id_cb, cutslist, vars_coeffs, vars_idx,
                                 nb_vars, sense, rhs)
end

function wbcm_attach_var_functor!(c_model::Ptr{Cvoid}, name::Symbol,
                                 c_fct::Ptr{Cvoid}, jump_var::Any)
  @bcm_ccall("attachVarFunctor", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cvoid}, Any),
                                        c_model, name, c_fct, jump_var)
end

function wbcm_attach_cstr_functor!(c_model::Ptr{Cvoid}, name::Symbol,
                                 c_fct::Ptr{Cvoid}, jump_cstr::Any)
  @bcm_ccall("attachCstrFunctor", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cvoid}, Any),
                                        c_model, name, c_fct, jump_cstr)
end

function wbcm_add_dyn_var!(c_model::Ptr{Cvoid}, col_id::Cint, name::Symbol,
                           var_bcid::Vector{Cint}, sp_bctype::Cint,
                           sp_bcid::Vector{Cint})
  @bcm_ccall("addDynVar", Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Cint,
                                 Ptr{Cint}),
                                 c_model, col_id, name, var_bcid, sp_bctype,
                                 sp_bcid)
end

function wbcm_add_dyn_cstr!(c_model::Ptr{Cvoid}, row_id::Cint, name::Symbol,
                            cstr_bcid::Vector{Cint})
  @bcm_ccall("addDynCstr", Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}),
                                  c_model, row_id, name, cstr_bcid)
end

function wbcm_get_cstr(c_model::Ptr{Cvoid}, row_id::Cint)
  @bcm_ccall("getCstr", Ptr{Cvoid}, (Ptr{Cvoid}, Cint), c_model, row_id)
end

function wbcs_add_cut!(cuts_list::Ptr{Cvoid}, cstr::Ptr{Cvoid})
  @bcs_ccall("addCut", Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cuts_list, cstr)
end

function wbcm_set_constr_bc_type!(cstr::Ptr{Cvoid}, bctype::Char)
  @bcm_ccall("setConstrBcType", Cint, (Ptr{Cvoid}, Cchar), cstr, bctype)
end

function wbcs_add_to_oracle_dual_sol!(c_model::Ptr{Cvoid}, c_solution::Ptr{Cvoid},
                                      cstr::Ptr{Cvoid})
  @bcs_ccall("addToOracleDualSol", Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
                                                c_model, c_solution, cstr)
end

function wbcm_add_cstr_term!(c_model::Ptr{Cvoid}, bc_cstr::Ptr{Cvoid},
                             col_id::Cint, coeff::Cdouble)
  @bcm_ccall("addCstrTerm", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cdouble),
                                   c_model, bc_cstr, col_id, coeff)
end

function wbcm_add_cstr_terms!(c_model::Ptr{Cvoid}, bc_cstr::Ptr{Cvoid}, 
                              cols_id::Vector{Cint}, coeff::Cdouble)
  @bcm_ccall("addCstrTerms", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cint}, Cint, Cdouble),
                c_model, bc_cstr, cols_id, length(cols_id), coeff)
end

function wbcm_add_membership_to_cstr!(c_model::Ptr{Cvoid}, memberships::Ptr{Cvoid},
                                     row_id::Cint, coeff::Cdouble)
  @bcm_ccall("addMembershipToCstr", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cdouble),
                                           c_model, memberships, row_id, coeff)
end

function wbcm_add_cstr_rhs!(bc_cstr::Ptr{Cvoid}, sense::Char, rhs::Cdouble)
  @bcm_ccall("addCstrRhs", Cvoid, (Ptr{Cvoid}, Cchar, Cdouble),
                                  bc_cstr, sense, rhs)
end
