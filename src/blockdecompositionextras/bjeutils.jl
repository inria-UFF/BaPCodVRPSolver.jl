isa_JuMPcstr(contnr) = isa(contnr, JuMP.ConstraintRef)
isa_JuMPvar(contnr) = isa(contnr, JuMP.Variable)


name(x::DynVariableDict) = name(x.jvar)

jumpmodel(x::DynVariable) = x.m
jumpmodel(c::DynConstraint) = c.m
jumpmodel(x::DynVariableDict) = jumpmodel(x.jvar)

function lookfornameofJuMPobj(m::JuMP.Model, c::JuMP.JuMPContainer)
  for (jm, jc) in m.objDict
    (jc == c) && return jm
  end
end
