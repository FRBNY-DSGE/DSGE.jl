abstract Scenario

function Base.show(io::IO, scen::Scenario)
    @printf io "Scenario: %s\n" scen.description
    @printf io "Targets: %s\n" scen.target_names
    @printf io "Instruments: %s\n" scen.instrument_names
end

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)

function get_scenario_input_file(m::AbstractModel, key::Symbol, vint::String)
    basename = lowercase(string(key)) * "_" * vint * ".jld"
    return inpath(m, "scenarios", basename)
end

function get_scenario_output_files(m::AbstractModel, key::Symbol, vint::String,
                                   output_vars::Vector{Symbol})
    filestring_addl = Vector{String}()
    push!(filestring_addl, "scen=" * string(key))
    push!(filestring_addl, "svin=" * vint)

    output_files = Dict{Symbol, String}()
    for var in output_vars
        basename = string(var) * ".jld"
        output_files[var] = workpath(m, "scenarios", basename, filestring_addl)
    end
    return output_files
end

function n_scenario_draws(m::AbstractModel, key::Symbol, vint::String)
    input_file = get_scenario_input_file(m, key, vint)
    draws = h5open(input_file, "r") do file
        dataset = HDF5.o_open(file, "arr")
        size(dataset)[3]
    end
    return draws
end

function load_scenario_targets!(scen::Scenario, path::String, draw_index::Int)
    raw_targets = h5read(path, "arr", (:, :, draw_index))
    target_inds = load(path, "target_indices")

    @assert collect(keys(target_inds)) == scen.target_names "Target indices in $path do not match target names in $(scen.key)"

    for (target_name, target_index) in target_inds
        scen.targets[target_name] = raw_targets[target_index, :]
    end

    return scen.targets
end

function targets_to_data(m::AbstractModel, scen::Scenario)
    df = DataFrame()
    horizons = n_target_horizons(scen)
    df[:date] = DSGE.quarter_range(date_forecast_start(m), DSGE.iterate_quarters(date_forecast_start(m), horizons - 1))
    for var in keys(m.observables)
        df[var] = if var in scen.target_names
            scen.targets[var]
        else
            fill(NaN, horizons)
        end
    end
    return df
end