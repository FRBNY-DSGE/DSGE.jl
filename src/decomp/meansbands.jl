function decomposition_means(m_new::M, m_old::M, input_type::Symbol,
                             cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol}; forecast_string_new = "", forecast_string_old = "",
                             verbose::Symbol = :low) where M<:AbstractModel
    # Print
    println(verbose, :low, )
    info_print(verbose, :low, "Computing means of forecast decomposition...")
    println(verbose, :low, "Start time: " * string(now()))
    begin_time = time_ns()

    input_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, classes, forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old)

    for class in classes
        print(verbose, :high, "Computing " * string(class) * "...")

        # Read metadata
        variable_names = JLD2.jldopen(input_files[Symbol(:decomptotal, class)]) do file
            class_long = get_class_longname(class)
            collect(keys(read(file, string(class_long) * "_indices")))
        end

        # Get to work!
        mapfcn = use_parallel_workers(m_new) ? pmap : map
        decomp_vec = mapfcn(var -> decomposition_means(m_new, m_old, input_type,
                                                       cond_new, cond_old, class, var, forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old),
                            variable_names)
        decomps = OrderedDict{Symbol, DataFrame}()
        for (var, decomp) in zip(variable_names, decomp_vec)
            decomps[var] = decomp
        end

        # Write to file
        output_file = get_decomp_mean_file(m_new, m_old, input_type, cond_new, cond_old, class, forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old)
        output_dir = dirname(output_file)
        isdir(output_dir) || mkpath(output_dir)
        JLD2.jldopen(output_file, "w") do file
            write(file, "decomps", decomps)
        end
        println(verbose, :high, "wrote " * basename(output_file))
    end

    # Print
    total_mb_time     = (time_ns() - begin_time)/1e9
    total_mb_time_min = total_mb_time/60

    println(verbose, :low, "\nTotal time to compute decomposition means and bands: " * string(total_mb_time_min) * " minutes")
    println(verbose, :low, "Computation of means and bands complete: " * string(now()))
end

function decomposition_means(m_new::M, m_old::M, input_type::Symbol,
                             cond_new::Symbol, cond_old::Symbol,
                             class::Symbol, var::Symbol; forecast_string_new = "", forecast_string_old = "") where M<:AbstractModel
    # Read in dates
    input_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, [class], forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old)
    input_file = input_files[Symbol(:decomptotal, class)]
    dates = JLD2.jldopen(input_file, "r") do file
        sort(collect(keys(read(file, "date_indices"))))
    end

    decomp = DataFrame(date = dates)

    for comp in [:data, :news, :shockdec, :dettrend, :para, :total]
        product = Symbol(:decomp, comp)

        input_file = input_files[Symbol(product, class)]
        JLD2.jldopen(input_file, "r") do file
            # Parse transform
            class_long = get_class_longname(class)
            transforms = read(file, string(class_long) * "_revtransforms")
            transform = parse_transform(transforms[var])

            # If shockdec, loop over shocks
            loopkeys = if comp == :shockdec
                shock_indices = read(file, "shock_indices")
                collect(keys(shock_indices))
            else
                [comp]
            end

            for key in loopkeys
                # Read in raw output: ndraws x nperiods
                decomp_series = if comp == :shockdec
                    read_forecast_series(file, class, product, var, key)
                else
                    read_forecast_series(file, class, product, var)
                end

                # Reverse transform
                transformed_decomp = scenario_mb_reverse_transform(decomp_series, transform, :forecast)

                # Compute mean and add to DataFrame
                decomp[key] = vec(mean(transformed_decomp, dims = 1))
            end
        end
    end

    return decomp
end
