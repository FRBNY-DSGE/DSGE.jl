function scenario_means_bands(m::AbstractModel, scenario_key::Symbol, scenario_vint::String,
                              output_vars::Vector{Symbol} = [:forecastutobs, :forecastutpseudo,
                                                             :forecastobs, :forecastpseudo,
                                                             :forecast4qobs, :forecast4qpseudo];
                              verbose::Symbol = :low,
                              kwargs...)
    # Print
    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        println()
        info("Computing means and bands for scenario = $scenario_key...")
        println("Start time: " * string(now()))
        println("Means and bands will be saved in " * workpath(m, "scenarios"))
        tic()
    end

    for output_var in output_vars
        if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
            print("Computing " * string(output_var) * "...")
        end

        # Compute means and bands
        input_file = get_scenario_mb_input_file(m, scenario_key, scenario_vint, output_var)
        mb = means_bands(scenario_key, scenario_vint, output_var, input_file)

        # Write to file
        output_file = get_scenario_mb_output_file(m, scenario_key, scenario_vint, output_var)
        output_dir = dirname(output_file)
        isdir(output_dir) || mkpath(output_dir)
        jldopen(output_file, "w") do file
            write(file, "mb", mb)
        end

        if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
            println("wrote " * basename(output_file))
        end
    end

    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        total_mb_time     = toq()
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute scenario means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function DSGE.means_bands(scenario_key::Symbol, scenario_vint::String, output_var::Symbol,
                          input_file::String; kwargs...)
    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read metadata
    metadata, mb_metadata =
        get_scenario_mb_metadata(scenario_key, scenario_vint, output_var, input_file)
    variable_names = collect(keys(mb_metadata[:indices]))
    date_list = collect(keys(mb_metadata[:date_inds]))

    # Get to work!
    mapfcn = use_parallel_workers(m) ? pmap : map
    mb_vec = pmap(var_name -> compute_scenario_means_bands(class, product, var_name, input_file),
                  variable_names)

    # Re-assemble pmap outputs
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()
    for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
        means[var_name] = var_means
        bands[var_name] = var_bands
        bands[var_name][:date] = date_list
    end

    return MeansBands(mb_metadata, means, bands; kwargs...)
end

function compute_scenario_means_bands(class::Symbol, product::Symbol, var_name::Symbol,
                                      filename::String;
                                      minimize::Bool = false,
                                      density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9])

    @assert product in [:forecast, :forecast4q, :forecastut] "Product can only be forecast(ut|4q)"

    # Read in everything from raw forecast file
    fcast_series, transform, var_ind, date_list = jldopen(filename, "r") do file
        # Read forecast output
        fcast_series = read_forecast_output(file, class, product, var_name)

        # Parse transform
        class_long = DSGE.get_class_longname(class)
        transforms = read(file, string(class_long) * "_revtransforms")
        transform = DSGE.parse_transform(transforms[var_name])

        # Get variable index
        indices = read(file, string(class_long) * "_indices")
        var_ind = indices[var_name]

        # Read date list
        date_list = collect(keys(read(file, "date_indices")))

        fcast_series, transform, var_ind, date_list
    end

    # Reverse transform
    if product == :forecast4q
        transform4q_scen = DSGE.get_scenario_transform4q(transform)

        y0s = if transform4q_scen == DSGE.loggrowthtopct_4q_approx
            # Sum growth rates y_{t-3}, y_{t-2}, y_{t-1}, and y_t
            zeros(3)
        elseif transform4q_scen == DSGE.logleveltopct_4q_approx
            # Divide log levels y_t by y_{t-4}
            zeros(4)
        else
            Float64[]
        end

        transformed_series = reverse_transform(fcast_series, transform4q_scen;
                                               fourquarter = true, y0s = y0s)
    elseif product == :forecast
        transform_scen = DSGE.get_scenario_transform(transform)
        transformed_series = reverse_transform(fcast_series, transform_scen; y0 = 0.0)

    else
        transformed_series = fcast_series
    end

    # Compute means and bands of transformed series
    means = vec(mean(transformed_series, 1))
    bands = find_density_bands(transformed_series, density_bands, minimize = minimize)
    return means, bands
end