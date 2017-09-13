function get_scenario_input_file(m::AbstractModel, key::Symbol, vint::String)
    basename = lowercase(string(key)) * "_" * vint * ".jld"
    return inpath(m, "scenarios", basename)
end

function get_scenario_filename(m::AbstractModel, key::Symbol,
                               vint::String, output_var::Symbol;
                               pathfcn::Function = rawpath,
                               fileformat = :jld)

    filestring_addl = Vector{String}()
    push!(filestring_addl, "scen=" * lowercase(string(key)))
    push!(filestring_addl, "svin=" * vint)

    basename = string(output_var) * "." * string(fileformat)
    return pathfcn(m, "scenarios", basename, filestring_addl)
end

function get_scenario_output_files(m::AbstractModel, key::Symbol, vint::String,
                                   output_vars::Vector{Symbol})
    output_files = Dict{Symbol, String}()
    for var in output_vars
        output_files[var] = get_scenario_filename(m, key, vint, var)
    end
    return output_files
end

function n_scenario_draws(m::AbstractModel, key::Symbol, vint::String)
    input_file = get_scenario_input_file(m, key, vint)
    draws = h5open(input_file, "r") do file
        dataset = HDF5.o_open(file, "arr")
        size(dataset)[1]
    end
    return draws
end

function load_scenario_targets!(scen::Scenario, path::String, draw_index::Int)
    raw_targets = squeeze(h5read(path, "arr", (draw_index, :, :)), 1)
    target_inds = load(path, "target_indices")

    @assert collect(keys(target_inds)) == scen.target_names "Target indices in $path do not match target names in $(scen.key)"

    for (target_name, target_index) in target_inds
        scen.targets[target_name] = raw_targets[target_index, :]
    end

    return scen.targets
end
