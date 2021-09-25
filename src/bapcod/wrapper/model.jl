function c_delete!(modelptr::Ptr{Cvoid})
    @bcm_ccall("delete", Cvoid, (Ptr{Cvoid},), modelptr)
end

function c_delete_solution!(sol)
  @bcsol_ccall("delete", Cvoid, (Ptr{Cvoid},), sol)
end

function init_model!(modelptr::Ptr{Cvoid}, nbrows::Cint, nbcols::Cint)
    @bcm_ccall("initModel", Cvoid, (Ptr{Cvoid}, Cint, Cint),
                          modelptr, nbrows, nbcols)
end

function new!()
    @bcsol_ccall("new", Ptr{Cvoid}, ())
end

function new!(param_file::String, print_param::Bool, int_obj::Bool, int_valued_bound::Bool, argc::Cint, argv::Array{String})
    @bcm_ccall("new", Ptr{Cvoid}, (Ptr{UInt8}, UInt8, UInt8, UInt8, Cint, Ptr{Ptr{UInt8}}),
                                param_file, print_param, int_obj, int_valued_bound, argc, argv)
end

function register_sub_problem!(mptr::Ptr{Cvoid}, subproblemtype::Cint, spmid::Array)
    @bcm_ccall("registerSubProblem", Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}),
                                        mptr, subproblemtype, spmid)
end

function set_art_cost_value!(mptr::Ptr{Cvoid}, acv::Cdouble)
    @bcm_ccall("setArtCostValue", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, acv)
end

function set_obj_lb!(mptr::Ptr{Cvoid}, lb::Cdouble)
    @bcm_ccall("setObjLb", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, lb)
end

function set_obj_magnitude!(mptr::Ptr{Cvoid}, ma::Cdouble)
    @bcm_ccall("setObjMagnitude", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, ma)
end

function set_obj_ub!(mptr::Ptr{Cvoid}, ub::Cdouble)
    @bcm_ccall("setObjUb", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, ub)
end

function wbcm_add_aggrsubprob_var_branching(c_model::Ptr{Cvoid}, varname::Symbol, highest_prio::Cdouble, prio::Cdouble, indicesignored::Cint, preproc::Bool)
  status = @bcm_ccall("addAggrSubProbVarBranching", Cint,
                        (Ptr{Cvoid}, Ptr{Cchar}, Cdouble, Cdouble, Cint, UInt8),
                        c_model, varname, highest_prio, prio, indicesignored, preproc)
  (status != 1) && error("Cannot add the aggrsubprob_var_branching.")
end

function wbcm_register_branching_expression(c_model::Ptr{Cvoid}, varname::Symbol,
        priority::Cdouble)
    arrayid = @bcm_ccall("registerBranchingExpression", Cint, 
                        (Ptr{Cvoid}, Ptr{Cchar}, Cdouble),
                        c_model, varname, priority)
    return arrayid
end

function wbcm_add_branching_expression(c_model::Ptr{Cvoid}, arrayid::Cint, 
        expmid::Vector{Cint}, varids::Vector{Cint}, coeffs::Vector{Cdouble}, l::Cint)
    status = @bcm_ccall("addBranchingExpression", Cint,
            (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
             c_model, arrayid, expmid, varids, coeffs, l)
    (status != 1) && error("Cannot add branching expression.")
end
