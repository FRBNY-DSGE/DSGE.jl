"""
```
get_scenario_input_file(m, scen::Scenario)
```

Get file name of raw scenario targets from `inpath(m, \"scenarios\")`.
"""
function get_scenario_input_file(m::AbstractModel, scen::Scenario)
    basename = string(scen.key) * "_" * scen.vintage * ".jld"
    return inpath(m, "scenarios", basename)
end

"""
```
n_scenario_draws(m, scen::Scenario)
```

Return the number of draws for `scen`, determined using
`get_scenario_input_file(m, scen)`.
"""
function n_scenario_draws(m::AbstractModel, scen::Scenario)
    input_file = get_scenario_input_file(m, scen)
    draws = h5open(input_file, "r") do file
        dataset = HDF5.o_open(file, "arr")
        size(dataset)[1]
    end
    return draws
end

"""
```
load_scenario_targets!(m, scen::Scenario, draw_index)
```

Add the targets from the `draw_index`th draw of the raw scenario targets to
`scen.targets`.
"""
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

"""
```
get_scenario_filename(m, scen::AbstractScenario, output_var;
    pathfcn = rawpath, fileformat = :jld, directory = "")
```

Get scenario file name of the form
`pathfcn(m, \"scenarios\", output_var * filestring * string(fileformat))`. If
`directory` is provided (nonempty), then the same file name in that directory
will be returned instead.
"""
function get_scenario_filename(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                               pathfcn::Function = rawpath,
                               fileformat::Symbol = :jld,
                               directory::String = "")
    filestring_addl = Vector{String}()
    if isa(scen, SingleScenario)
        push!(filestring_addl, "scen=" * string(scen.key))
    elseif isa(scen, ScenarioAggregate)
        push!(filestring_addl, "sagg=" * string(scen.key))
    end
    push!(filestring_addl, "svin=" * scen.vintage)

    base = string(output_var) * "." * string(fileformat)
    path = pathfcn(m, "scenarios", base, filestring_addl)
    if !isempty(directory)
        path = joinpath(directory, basename(path))
    end
    return path
end

"""
```
get_scenario_output_files(m, scen::SingleScenario, output_vars)
```

Return a `Dict{Symbol, String}` mapping `output_vars` to the raw simulated
scenario outputs for `scen`.
"""
function get_scenario_output_files(m::AbstractModel, scen::SingleScenario,
                                   output_vars::Vector{Symbol})
    output_files = Dict{Symbol, String}()
    for var in output_vars
        output_files[var] = get_scenario_filename(m, scen, var)
    end
    return output_files
end

"""
```
get_scenario_mb_input_file(m, scen::AbstractScenario, output_var)
```

Call `get_scenario_filename` while replacing `forecastut` and `forecast4q` in
`output_var` with `forecast`.
"""
function get_scenario_mb_input_file(m::AbstractModel, scen::AbstractScenario, output_var::Symbol)
    input_file = get_scenario_filename(m, scen, output_var)
    input_file = replace(input_file, "forecastut", "forecast")
    input_file = replace(input_file, "forecast4q", "forecast")
    return input_file
end

"""
```
get_scenario_mb_output_file(m, scen::AbstractScenario, output_var;
    directory = "")
```

Call `get_scenario_filename` while tacking on `\"mb\"` to the front of the base
file name.
"""
function get_scenario_mb_output_file(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                                     directory::String = "")
    fullfile = get_scenario_filename(m, scen, output_var, pathfcn = workpath, directory = directory)
    joinpath(dirname(fullfile), "mb" * basename(fullfile))
end

"""
```
get_scenario_mb_metadata(m, scen::SingleScenario, output_var)

get_scenario_mb_metadata(m, agg::ScenarioAggregate, output_var)
```

Return the `MeansBands` metadata dictionary for `scen`.
"""
function get_scenario_mb_metadata(m::AbstractModel, scen::SingleScenario, output_var::Symbol)
    forecast_output_file = get_scenario_mb_input_file(m, scen, output_var)
    _, mb_metadata = DSGE.get_mb_metadata(:mode, :none, output_var, forecast_output_file)
    mb_metadata[:scenario_key] = scen.key
    mb_metadata[:scenario_vint] = scen.vintage

    return mb_metadata
end

function get_scenario_mb_metadata(m::AbstractModel, agg::ScenarioAggregate, output_var::Symbol)
    forecast_output_file = get_scenario_mb_input_file(m, agg.scenario_groups[1][1], output_var)
    _, metadata = DSGE.get_mb_metadata(:mode, :none, output_var, forecast_output_file)

    # Initialize start and end date
    start_date = date_forecast_start(m)
    end_date   = maximum(keys(metadata[:date_inds]))

    for scen in vcat(agg.scenario_groups...)
        forecast_output_file = get_scenario_mb_input_file(m, scen, output_var)
        scen_dates = jldopen(forecast_output_file, "r") do file
            read(file, "date_indices")
        end

        # Throw error if start date for this scenario doesn't match
        if map(reverse, scen_dates)[1] != start_date
            error("All scenarios in agg must start from the same date")
        end

        # Update end date if necessary
        end_date = max(end_date, maximum(keys(metadata[:date_inds])))
    end

    dates = DSGE.quarter_range(start_date, end_date)
    metadata[:date_inds] = OrderedDict{Date, Int}(d => i for (i, d) in enumerate(dates))
    metadata[:scenario_key] = agg.key
    metadata[:scenario_vint] = agg.vintage

    return metadata
end

"""
```
read_scenario_output(m, scen::SingleScenario, class, product, var_name)

read_scenario_output(m, agg::ScenarioAggregate, class, product, var_name)
```

Given either `scen` or `agg`, read in and return all draws of and the
appropriate reverse transform for `var_name`.
"""
function read_scenario_output(m::AbstractModel, scen::SingleScenario, class::Symbol, product::Symbol,
                              var_name::Symbol)
    # Get filename
    filename = get_scenario_mb_input_file(m, scen, Symbol(product, class))

    fcast_series, transform = jldopen(filename, "r") do file
        # Read forecast outputs
        fcast_series = read_forecast_output(file, class, product, var_name)

        # Parse transform
        class_long = DSGE.get_class_longname(class)
        transforms = read(file, string(class_long) * "_revtransforms")
        transform = DSGE.parse_transform(transforms[var_name])

        fcast_series, transform
    end
end

function read_scenario_output(m::AbstractModel, agg::ScenarioAggregate, class::Symbol,
                              product::Symbol, var_name::Symbol)
    # Aggregate scenarios
    ngroups = length(agg.scenario_groups)
    agg_draws = Vector{Matrix{Float64}}(ngroups)

    for (i, (group, pct)) in enumerate(zip(agg.scenario_groups, agg.proportions))
        # Read in draws for each scenario in group
        nscenarios = length(group)
        group_draws = Vector{Matrix{Float64}}(nscenarios)
        for (j, scen) in enumerate(group)
            filename = get_scenario_mb_input_file(m, scen, Symbol(product, class))
            group_draws[j] = jldopen(filename, "r") do file
                read_forecast_output(file, class, product, var_name)
            end
        end

        # Stack all scenario draws in this group together
        all_group_draws = cat(1, group_draws...)

        # Sample
        actual_ndraws = size(all_group_draws, 1)
        desired_ndraws = convert(Int, round(pct * agg.total_draws))

        sampled_inds = if agg.replace
            sample(1:actual_ndraws, desired_ndraws, replace = true)
        else
            if desired_ndraws == 0
                Int[]
            else
                quotient  = convert(Int, floor(actual_ndraws / desired_ndraws))
                remainder = actual_ndraws % desired_ndraws
                vcat(repmat(1:actual_ndraws, quotient),
                     sample(1:actual_ndraws, remainder, replace = false))
            end
        end
        sort!(sampled_inds)
        agg_draws[i] = all_group_draws[sampled_inds, :]
    end

    fcast_series = cat(1, agg_draws...)

    # Parse transform
    filename = get_scenario_mb_input_file(m, agg.scenario_groups[1][1], Symbol(product, class))
    transform = jldopen(filename, "r") do file
        class_long = DSGE.get_class_longname(class)
        transforms = read(file, string(class_long) * "_revtransforms")
        DSGE.parse_transform(transforms[var_name])
    end

    return fcast_series, transform
end

"""
```
read_scenario_mb(m, scen::AbstractScenario, output_var; directory = "")
```

Read in an alternative scenario `MeansBands` object.
"""
function read_scenario_mb(m::AbstractModel, scen::AbstractScenario, output_var::Symbol;
                          directory::String = "")
    filepath = get_scenario_mb_output_file(m, scen, output_var, directory = directory)
    read_mb(filepath)
end