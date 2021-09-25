module BlockSolverExtraInterface

adddynamicvariable!() = nothing
export adddynamicvariable!

adddynamicconstraint!() = nothing
export adddynamicconstraint!

adddynamiccut!() = nothing
export adddynamiccut!

set_var_generic!() = nothing
export set_var_generic!

set_cstr_generic!() = nothing
export set_cstr_generic!

set_cstr_preproc!() = nothing
export set_cstr_preproc!

###### RCSP ####

set_rcsp_oracles!() = nothing
export set_rcsp_oracles!

set_rcsp_cuts!() = nothing
export set_rcsp_cuts!

set_rcsp_branching!() = nothing
export set_rcsp_branching!

send_rcsp_to_solver!() = nothing
export send_rcsp_to_solver!

##############
addtermtodynamiccstr!() = nothing
export addtermtodynamiccstr!

addtermstodynamiccstr!() = nothing
export addtermstodynamiccstr!

addrhstodynamiccstr!() = nothing
export addrhstodynamiccstr!

addentrytodynamiccol!() = nothing
export addentrytodynamiccol!

setobjvaluetodynamiccol!() = nothing
export setobjvaluetodynamiccol!

getcurcostofdynvar() = nothing
export getcurcostofdynvar

getdualofdyncstr() = nothing
export getdualofdyncstr

#########
getvalueofdynvar() = nothing
export getvalueofdynvar end
