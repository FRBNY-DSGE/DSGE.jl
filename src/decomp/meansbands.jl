function decomposition_means(m_new::M, m_old::M, input_type::Symbol,
                             cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol};
                             verbose::Symbol = :low) where M<:AbstractModel
    # Print
    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        println()
        info("Computing means of forecast decomposition...")
        println("Start time: " * string(now()))
        tic()
    end

    input_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, classes)

    for class in classes
        if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
            print("Computing " * string(class) * "...")
        end

        # Read metadata
        variable_names = jldopen(input_files[Symbol(:decomptotal, class)]) do file
            class_long = DSGE.get_class_longname(class)
            collect(keys(read(file, string(class_long) * "_indices")))
        end

        # Get to work!
        mapfcn = use_parallel_workers(m_new) ? pmap : map
        decomp_vec = mapfcn(var -> decomposition_means(m_new, m_old, input_type,
                                                       cond_new, cond_old, class, var),
                            variable_names)
        decomps = OrderedDict{Symbol, DataFrame}()
        for (var, decomp) in zip(variable_names, decomp_vec)
            decomps[var] = decomp
        end

        # Write to file
        output_file = get_decomp_mean_file(m_new, m_old, input_type, cond_new, cond_old, class)
        output_dir = dirname(output_file)
        isdir(output_dir) || mkpath(output_dir)
        jldopen(output_file, "w") do file
            write(file, "decomps", decomps)
        end

        if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
            output_file = get_scenario_mb_output_file(m, scen, output_var)
            println("wrote " * basename(output_file))
        end
    end

    # Print
    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        total_mb_time     = toq()
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute decomposition means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function decomposition_means(m_new::M, m_old::M, input_type::Symbol,
                             cond_new::Symbol, cond_old::Symbol,
                             class::Symbol, var::Symbol) where M<:AbstractModel
    # Read in dates
    input_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, [class])
    input_file = input_files[Symbol(:decomptotal, class)]
    dates = jldopen(input_file, "r") do file
        sort(collect(keys(read(file, "date_indices"))))
    end

    decomp = DataFrame(date = dates)

    for comp in [:data, :news, :shockdec, :dettrend, :para, :total]
        product = Symbol(:decomp, comp)

        input_file = input_files[Symbol(product, class)]
        jldopen(input_file, "r") do file
            # Parse transform
            class_long = DSGE.get_class_longname(class)
            transforms = read(file, string(class_long) * "_revtransforms")
            transform = DSGE.parse_transform(transforms[var])

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
                    DSGE.read_forecast_series(file, class, product, var, key)
                else
                    DSGE.read_forecast_series(file, class, product, var)
                end

                # Reverse transform
                transformed_decomp = DSGE.scenario_mb_reverse_transform(decomp_series, transform, :forecast)

                # Compute mean and add to DataFrame
                decomp[key] = vec(mean(transformed_decomp, 1))
            end
        end
    end

    return decomp
end
