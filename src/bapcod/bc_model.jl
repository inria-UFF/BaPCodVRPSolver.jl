mutable struct BcModel
  ptr::Ptr{Cvoid} # BaPCod problem
  nb_vars::Int
  nb_cstrs::Int
  solution::Ptr{Cvoid} # Solution to the problem

  function BcModel(;param_file = "config/bcParameters.cfg",
                    print_param = true, integer_objective = false, 
                    baptreedot_file = "BaPTree.dot", user_params = "")
    # TODO check
    (!isfile(param_file)) &&
      error("Parameter file $param_file does not exist. Create a config file or
      copy one from BaPCod.jl/configexamples/")
    user_params_array = map(String, split(user_params, " ")) # to create an array of Strings instead of SubStrings
    argv = ["empty", "-t", baptreedot_file]
    if length(user_params_array) > 0
      argv = vcat(argv, user_params_array)
    end
    argc = Cint(length(argv))
    modelptr = new!(param_file, print_param, integer_objective, integer_objective, argc, argv)
    solptr = new!()
    model = new(modelptr, 0, 0, solptr)
    finalizer(delete_model, model)
    return model
  end
end

function delete_model(model::BcModel)
  if model.ptr != C_NULL
    c_delete!(model.ptr)
    model.ptr = C_NULL
  end

  if model.solution != C_NULL
    c_delete_solution!(model.solution)
    model.solution = C_NULL
  end
  return
end

function c_init_model!(model::BcModel, nbrows, nbcols)
  init_model!(model.ptr, Cint(nbrows), Cint(nbcols))
  model.nb_vars = nbcols
  model.nb_cstrs = nbrows

  c_set_artcostvalue(model, 2.3)
end

nb_vars(model::BcModel) = model.nb_vars
nb_cstrs(model::BcModel) = model.nb_vars

function c_set_obj_sense!(model::BcModel, sense)
  (sens == :Max) && error("BaPCod framework does not support maximisation.")
end

# BaPCod framework does not support maximisation.
c_get_obj_sense(model::BcModel) = :Min

c_set_objlb(m::BcModel, lb) =
  set_obj_lb!(m.ptr, Cdouble(lb))

c_set_objub(m::BcModel, ub) =
  set_obj_ub!(m.ptr, Cdouble(ub))

c_set_objmagnitude(m::BcModel, ma) =
  set_obj_magnitude!(m.ptr, Cdouble(ma))

c_set_artcostvalue(m::BcModel, acv) =
  set_art_cost_value!(m.ptr, Cdouble(acv))

function c_register_subproblems(m::BcModel, spids)
  for (spid, sptype) in spids
    spmid = Array{Cint}(undef, 8)
    from_index_to_BaPCodindex(spid, spmid)
    subproblemtype = Cint(sptype_to_int(sptype))
    register_sub_problem!(m.ptr, subproblemtype, spmid)
  end
end
