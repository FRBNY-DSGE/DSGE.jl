"""
```
scenario_means_bands(m, scen::AbstractScenario,
    output_vars = [:forecastutobs, :forecastobs, :forecast4qobs,
                   :forecastutpseudo, :forecastpseudo, :forecast4qpseudo];
    verbose = :low, kwargs...)

scenario_means_bands(m, scen::AbstractScenario, output_var; kwargs...)

scenario_means_bands(m, scen::AbstractScenario, output_var, var_name; kwargs...)
```

Compute means and bands for model `m` and scenario `scen`. Keyword arguments are
the same as for `compute_scenario_means_bands`.
"""
function scenario_means_bands(m::AbstractModel, scen::AbstractScenario,
                              output_vars::Vector{Symbol} = [:forecastutobs, :forecastutpseudo,
                                                             :forecastobs, :forecastpseudo,
                                                             :forecast4qobs, :forecast4qpseudo];
                              verbose::Symbol = :low,
                              kwargs...)
    # Print
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        info("Computing means and bands for scenario = " * string(scen.key) * "...")
        println("Start time: " * string(now()))
        println("Means and bands will be saved in " * workpath(m, "scenarios"))
        tic()
    end

    for output_var in output_vars
        if VERBOSITY[verbose] >= VERBOSITY[:high]
            print("Computing " * string(output_var) * "...")
        end

        # Compute means and bands
        mb = scenario_means_bands(m, scen, output_var; kwargs...)

        # Write to file
        output_file = get_scenario_mb_output_file(m, scen, output_var)
        output_dir = dirname(output_file)
        isdir(output_dir) || mkpath(output_dir)
        jldopen(output_file, "w") do file
            write(file, "mb", mb)
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("wrote " * basename(output_file))
        end
    end

    # Print
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = toq()
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute scenario means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function scenario_means_bands(m::AbstractModel, scen::AbstractScenario, output_var::Symbol; kwargs...)
    # Read metadata
    metadata = get_scenario_mb_metadata(m, scen, output_var)
    date_list = collect(keys(metadata[:date_inds]))
    variable_names = collect(keys(metadata[:indices]))

    # Get to work!
    mapfcn = use_parallel_workers(m) ? pmap : map
    mb_vec = pmap(var_name -> scenario_means_bands(m, scen, output_var, var_name),
                  variable_names)

    # Re-assemble pmap outputs
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()
    for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
        means[var_name] = var_means
        bands[var_name] = var_bands
        bands[var_name][:date] = date_list
    end

    return MeansBands(metadata, means, bands; kwargs...)
end

function scenario_means_bands(m::AbstractModel, scen::AbstractScenario, output_var::Symbol,
                              var_name::Symbol; kwargs...)
    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in scenario draws
    fcast_series, transform = read_scenario_output(m, scen, class, product, var_name)

    # Compute means and bands
    compute_scenario_means_bands(fcast_series, transform, product; kwargs...)
end

"""
```
compute_scenario_means_bands(fcast_series, transform, product;
    minimize = false, density_bands = [0.5, 0.6, 0.7, 0.8, 0.9])
```

### Inputs

- `fcast_series::Matrix{Float64}`: `ndraws` x `horizon` matrix of untransformed
  scenario forecasts
- `transform::Function`: reverse transform (possibly no transform at all or 4Q)
- `product::Symbol`: must be one of `:forecast`, `:forecastut`, or `:forecast4q`

### Keyword Arguments

- `minimize::Bool`: if `true`, choose shortest interval, otherwise just chop off
  lowest and highest (percent/2)
- `density_bands::Vector{Float64}`: a vector of percent values (between 0 and 1) for
  which to compute density bands
"""
function compute_scenario_means_bands(fcast_series::Matrix{Float64}, transform::Function,
                                      product::Symbol; minimize::Bool = false,
                                      density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9])

    @assert product in [:forecast, :forecastut, :forecast4q] "Product can only be forecast(ut|4q)"

    # Reverse transform
    if product == :forecast4q
        transform4q_scen = get_scenario_transform4q(transform)

        y0s = if transform4q_scen == loggrowthtopct_4q_approx
            # Sum growth rates y_{t-3}, y_{t-2}, y_{t-1}, and y_t
            zeros(3)
        elseif transform4q_scen == logleveltopct_4q_approx
            # Divide log levels y_t by y_{t-4}
            zeros(4)
        else
            Float64[]
        end

        transformed_series = reverse_transform(fcast_series, transform4q_scen;
                                               fourquarter = true, y0s = y0s)
    elseif product == :forecast
        transform_scen = get_scenario_transform(transform)
        transformed_series = reverse_transform(fcast_series, transform_scen; y0 = 0.0)

    else
        transformed_series = fcast_series
    end

    # Compute means and bands of transformed series
    means = vec(mean(transformed_series, 1))
    bands = find_density_bands(transformed_series, density_bands, minimize = minimize)
    return means, bands
end