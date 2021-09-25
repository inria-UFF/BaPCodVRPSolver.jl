# Entry method for constraints registration with cstr names and indexes
function c_register_cstrs(m::BcModel, A, lb, ub, cstrs_decomposition)
  cstr_bcid = Array{Cint}(undef, 8)
  sp_bcid = Array{Cint}(undef, 8)

  for (row_id, (name, c_id, sp_type, sp_id)) in enumerate(cstrs_decomposition)
    # BaPCod needs an idnex
    c_id = (c_id == nothing) ? (0,) : c_id
    sp_id = (sp_id == nothing) ? (0,) : sp_id

    from_index_to_BaPCodindex(c_id, cstr_bcid)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    # Register the constraint
    c_register_constraint(m.ptr, name, Cint(row_id - 1), sp_type, sp_bcid, cstr_bcid)
  end
  c_init_constraints(m.ptr, A, lb, ub)
end

# Entry method for cstrs registration without cstrs names and indexes
function c_register_cstrs(m::BcModel, A, lb, ub)
  nb_rows, nb_cols = size(A)
  empty = [0, -1, -1, -1, -1, -1, -1, -1]
  for row_id in 1:nb_rows
    c_register_constraint(m.ptr, :cstr, Cint(row_id - 1), :MIP, empty, empty)
  end
  c_init_constraints(m.ptr, A, lb, ub)
end

# Registering constraint in BaPCod
function c_register_constraint(mptr::Ptr{Cvoid}, name::Symbol, row_id::Cint,
                    sp_type::Symbol, sp_mid::Vector{Int}, cstr_mid::Vector{Int})
  sp_mid2 = [Cint(s) for s in sp_mid]
  cstr_mid2 = [Cint(s) for s in cstr_mid]
  c_register_constraint(mptr, name, row_id, sp_type, sp_mid2, cstr_mid2)
end

function c_register_constraint(mptr::Ptr{Cvoid}, name::Symbol, row_id::Cint,
                    sp_type::Symbol, sp_mid::Vector{Cint}, cstr_mid::Vector{Cint})
  if sp_type != :ALL
    stat = register_cstr!(mptr, name, row_id, sptype_to_int(sp_type), sp_mid, cstr_mid)
    (stat != 1) && error("Cannot init constraint $name at row $row_id.")
  else
    stat = register_generic_cstr!(mptr, name, row_id)
    (stat != 1) && error("Cannot init the generic constraint $name at row $row_id")
  end
end

function c_init_constraints(mptr::Ptr{Cvoid}, A, lb, ub)
  coeffs = Vector{Cint}(A.colptr .- 1)
  rows_id = Vector{Cint}(rowvals(A) .- 1)
  stat = init_cstrs!(mptr, coeffs, rows_id, nonzeros(A), lb, ub)
  (stat != 1) && error("Cannot create constraints.")
end

function c_cstrs_used_in_preproc(m::BcModel, cstrs_used_in_preproc)
  sp_bcid = Array{Cint}(undef, 8)
  for ((cstr_name, (sp_type, sp_id)), used) in cstrs_used_in_preproc
    sp_bctype = sptype_to_int(sp_type)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    status = c_cstr_used_in_preprocessing!(m.ptr, cstr_name, sp_bctype, sp_bcid, used)
    (status != 1) && error("Cannot ")
  end
end
