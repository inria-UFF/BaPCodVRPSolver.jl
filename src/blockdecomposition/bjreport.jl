mutable struct VarCstrReport
  cstrs_report
  vars_report
end
VarCstrReport() = VarCstrReport(nothing, nothing)

contains(p, s) = findnext(s, p, 1) != nothing

function report_cstrs_and_vars!(m::JuMP.Model)
  nbrows = length(m.linconstr)
  nbcols = length(m.colNames)
  m.ext[:varcstr_report].cstrs_report = Array{Tuple{Symbol, Union{Tuple,Cvoid}}}(undef, nbrows)
  m.ext[:varcstr_report].vars_report = Array{Tuple{Symbol, Union{Tuple,Cvoid}}}(undef, nbcols)
  report_names_and_indexes!(m.ext[:varcstr_report], m.objDict)
  check_and_name_anonymous(m)
end

function name_anonymous(report, name)
  anonymous = 0
  for i in 1:length(report)
    if !isassigned(report, i)
      anonymous += 1
      report[i] = (name,(anonymous,))
    end
  end
end

function check_for_anonymous(report)
  for i in 1:length(report)
    if !isassigned(report, i)
      info = "Make sure that all variables and constraints have a name."
      errmsg = "BlockDecomposition does not support anonymous variables or constraints."
      bjerror(info, errmsg)
    end
  end
end

function check_and_name_anonymous(model)
  if use_DantzigWolfe(model)
    name_anonymous(model.ext[:varcstr_report].cstrs_report, :anonymous_cstr)
    check_for_anonymous(model.ext[:varcstr_report].vars_report)
  else
    check_for_anonymous(model.ext[:varcstr_report].cstrs_report)
    check_for_anonymous(model.ext[:varcstr_report].vars_report)
  end
end

function report_names_and_indexes!(report, dict)
  for collection in dict
    name = string(collection.first)
    # Is it a JuMP Container ? or an array ?
    isjumpcontnr = isa_jumpcontnr(collection.second)
    isarray = isa_array(collection.second)
    if isjumpcontnr || isarray
      contains_jumpcstr(collection.second) && add_names_and_indexes!(report.cstrs_report, collection)
      contains_jumpvar(collection.second) && add_names_and_indexes!(report.vars_report, collection)
    end
    # Is it a single constraint or a single variable ?
    isjumpcstr = isa_jumpcstr(collection.second)
    isjumpvar = isa_jumpvar(collection.second)
    isjumpcstr && add_name_and_index!(report.cstrs_report, collection)
    isjumpvar && add_name_and_index!(report.vars_report, collection)
    # error
    (!isjumpcontnr && !isarray && !isjumpcstr && !isjumpvar) && bjerror("Unsupported type : collection name = $name and type = $(typeof(collection.second))") #TODO
  end
end

function getpos(vc) :: Integer
  if isa_jumpcstr(vc)
    return vc.idx
  end
  if isa_jumpvar(vc)
    return vc.col
  end
end

function add_names_and_indexes!(report, collection)
  name = collection.first
  for index in keys(collection.second)
      pos = getpos(collection.second[index...])
      report[pos] = (name, index)
  end
end

function add_name_and_index!(report, singleton)
  name = singleton.first
  pos = getpos(singleton.second)
  report[pos] = (name, nothing)
end
