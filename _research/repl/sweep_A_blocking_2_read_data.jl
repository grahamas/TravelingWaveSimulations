(example, mdb) = load_data(datadir(), "reduced_line_dos_effectively_sigmoid")
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)
@show mods

is_range(x::Nothing) = false
is_range(x::AbstractRange) = true
is_range(x::Number) = false

all_mod_names = keys(mods) |> collect
all_mod_values = values(mods) |> collect
varied_mods = is_range.(all_mod_values)
fixed_mods = .!varied_mods
fixed_mods_dict = Dict(name => mods[name] for name in all_mod_names[fixed_mods])

mod_names = all_mod_names[varied_mods] |> sort
mod_values = [mods[name] for name in mod_names]
mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))

mod_array(T) = NamedAxisArray{Tuple(mod_names)}(Array{Union{Bool,Missing}}(undef, length.(mod_values)...), Tuple(mod_values))

first_result = MultiDBRowIter(mdb) |> first
classification_names = fieldnames(ExecutionClassifications)
classifications_A = NamedTuple{Tuple(classification_names)}([mod_array(Bool) for _ in classification_names]) 

for (this_mod, this_result) in MultiDBRowIter(mdb)
    exec_classification = this_result[:wave_properties]
    if exec_classification.has_propagation == true
        @show "HEI"
    end
    classifications_A_idx = this_mod[mod_names]
    if exec_classification === missing
        for name in classification_names
            classifications_A[name][classifications_A_idx...] = missing
        end
    else
        for name in classification_names
            classifications_A[name][classifications_A_idx...] = getproperty(exec_classification, name)
        end
    end
end
