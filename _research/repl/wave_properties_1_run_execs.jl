line_example = get_example("reduced_line_dos_effectively_sigmoid")

mods_no_prop_nt = (Aee=90.0,Aei=160.0, Aie=132.0, Aii=30.0)
(name_no_prop, exec_no_prop) = execute_single_modification(line_example, mods_no_prop_nt)

mods_front_nt = (Aee=167.0,Aei=95.0, Aie=132.0, Aii=30.0)
(name_front, exec_front) = execute_single_modification(line_example, mods_front_nt)
