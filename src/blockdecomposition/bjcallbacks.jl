function addfacultativecutcallback(model::JuMP.Model, f::Function)
  if model.ext[:facultative_cuts_cb] != nothing
    error("Facultative cuts callback already defined.")
  end
  model.ext[:facultative_cuts_cb] = f
  return
end

function addcorecutcallback(model::JuMP.Model, f::Function)
  if model.ext[:core_cuts_cb] != nothing
    error("Core cuts callback alread defined.")
  end
  model.ext[:core_cuts_cb] = f
  return
end

function getvalue(cb::MathProgBase.MathProgCallbackData, var::JuMP.Variable)
  colid = Cint(var.col)
  if !applicable(getvalvariable, cb, colid)
    error("Cannot get value of $var")
  end
  return getvalvariable(cb, colid)
end
