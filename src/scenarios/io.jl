function get_scenario_input_file(m::AbstractModel, scen::Scenario)
    basename = lowercase(string(scen.key)) * "_" * scen.vintage * ".jld"
    return inpath(m, "scenarios", basename)
end

function get_scenario_filename(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                               pathfcn::Function = rawpath,
                               fileformat::Symbol = :jld,
                               directory::String = "")
    filestring_addl = Vector{String}()
    if isa(scen, SingleScenario)
        push!(filestring_addl, "scen=" * lowercase(string(scen.key)))
    elseif isa(scen, ScenarioAggregate)
        push!(filestring_addl, "sagg=" * lowercase(string(scen.key)))
    end
    push!(filestring_addl, "svin=" * scen.vintage)

    base = string(output_var) * "." * string(fileformat)
    path = pathfcn(m, "scenarios", base, filestring_addl)
    if !isempty(directory)
        path = joinpath(directory, basename(path))
    end
    return path
end

function get_scenario_output_files(m::AbstractModel, scen::SingleScenario,
                                   output_vars::Vector{Symbol})
    output_files = Dict{Symbol, String}()
    for var in output_vars
        output_files[var] = get_scenario_filename(m, scen, var)
    end
    return output_files
end

function n_scenario_draws(m::AbstractModel, scen::SingleScenario)
    input_file = get_scenario_input_file(m, scen)
    draws = h5open(input_file, "r") do file
        dataset = HDF5.o_open(file, "arr")
        size(dataset)[1]
    end
    return draws
end

function load_scenario_targets!(m::AbstractModel, scen::Scenario, draw_index::Int)
    path = get_scenario_input_file(m, scen)
    raw_targets = squeeze(h5read(path, "arr", (draw_index, :, :)), 1)
    target_inds = load(path, "target_indices")

    @assert collect(keys(target_inds)) == scen.target_names "Target indices in $path do not match target names in $(scen.key)"

    for (target_name, target_index) in target_inds
        scen.targets[target_name] = raw_targets[target_index, :]
    end

    return scen.targets
end

function get_scenario_mb_input_file(m::AbstractModel, scen::AbstractScenario, output_var::Symbol)
    input_file = get_scenario_filename(m, key, vint, output_var)
    input_file = replace(input_file, "forecastut", "forecast")
    input_file = replace(input_file, "forecast4q", "forecast")
    return input_file
end

function get_scenario_mb_output_file(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                                     directory::String = "")
    fullfile = get_scenario_filename(m, scen, output_var, pathfcn = workpath, directory = directory)
    joinpath(dirname(fullfile), "mb" * basename(fullfile))
end

"""
```
get_scenario_mb_metadata(m, scen, output_var)
```

Returns the `MeansBands` metadata dictionary for the `SingleScenario` `scen`.
"""
function get_scenario_mb_metadata(m::AbstractModel, scen::SingleScenario, output_var::Symbol)
    forecast_output_file = get_scenario_mb_input_file(m, scen, output_var)
    _, mb_metadata = DSGE.get_mb_metadata(:mode, :none, output_var, forecast_output_file)
    mb_metadata[:scenario_key] = scen.key
    mb_metadata[:scenario_vint] = scen.vintage

    return mb_metadata
end

function get_scenario_mb_metadata(m::AbstractModel, agg::ScenarioAggregate, output_var::Symbol)
    _, metadata = get_scenario_mb_input_file(m, agg.scenario_groups[1][1], output_var)
    start_date = date_forecast_start(m)

    for scen in vcat(agg.scenario_groups...)
        forecast_output_file = get_scenario_mb_input_file(m, scen, output_var)
        scen_dates = jldopen(forecast_output_file, "r") do file
            read(file, "date_indices")
        end

        # Throw error if start date for this scenario doesn't match
        if map(reverse, scen_dates)[1] != start_date
            error("All scenarios in agg must start from the same date")
        end

        merge!(metadata[:date_inds], scen_dates)
    end

    metadata[:scenario_key] = agg.key
    metadata[:scenario_vint] = agg.vintage

    return metadata
end

"""
```
read_scenario_mb(m, key, vint, output_var; directory = "")
```

Read in an alternative scenario `MeansBands` object.
"""
function read_scenario_mb(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                          directory::String = "")
    filepath = get_scenario_mb_output_file(m, scen, output_var, directory = directory)
    read_mb(filepath)
end