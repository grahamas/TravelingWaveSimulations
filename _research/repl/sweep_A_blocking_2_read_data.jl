(example, mdb) = load_data(datadir(), "reduced_line_dos_effectively_sigmoid")
#example_name = get_example_name(mdb.fns[1])
#sim_name = get_sim_name(mdb.fns[1])

GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)

all_mod_names = keys(mods) |> collect
all_mod_values = values(mods) |> collect
varied_mods = length.(all_mod_values) .> 1
fixed_mods = length.(all_mod_values) .== 1
fixed_mods_dict = Dict(name => mods[name] for name in all_mod_names[fixed_mods])

mod_names = all_mod_names[varied_mods] |> sort
mod_values = [mods[name] for name in mod_names]
mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))

mod_array(T) = AxisArray(Array{Union{Bool,Missing}}(undef, length.(mod_values)...), Tuple(mod_values))

first_result = MultiDBRowIter(mdb) |> first
classification_names = fieldnames(ExecutionClassifications)
A = NamedTuple{Tuple(classification_names)}([mod_array(Bool) for _ in classification_names]) 

for (this_mod, this_result) in MultiDBRowIter(mdb)
    exec_classification = this_result[:wave_properties]
    A_idx = this_mod[mod_names]
    if exec_classification === missing
        for name in classification_names
            A[name][A_idx...] = missing
        end
    else
        for name in classification_names
            A[name][A_idx...] = getproperty(exec_classification, name)
        end
    end
end