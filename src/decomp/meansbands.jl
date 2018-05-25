function decomposition_means(m_new::AbstractModel, m_old::AbstractModel,
                             input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                             classes::Vector{Symbol}; verbose::Symbol = :low)
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
        decomps = DataStructures.OrderedDict{Symbol, DataFrame}()
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

        println("\nTotal time to compute scenario means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function decomposition_means(m_new::AbstractModel, m_old::AbstractModel,
                             input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                             class::Symbol, var::Symbol)
    # Read in dates
    input_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, [class])
    input_file = input_files[Symbol(:decomptotal, class)]
    dates = jldopen(input_file, "r") do file
        sort(collect(keys(read(file, "date_indices"))))
    end

    decomp = DataFrame(date = dates)

    # Non-individual shock components are ndraws x nperiods
    for comp in [:state, :shock, :data, :param, :total]
        product = Symbol(:decomp, comp)

        # Read in raw output
        input_file = input_files[Symbol(product, class)]
        decomp_series, transform = jldopen(input_file, "r") do file
            # Read in variable: ndraws x nperiods
            decomp_series = DSGE.read_forecast_series(file, class, product, var)

            # Parse transform
            class_long = DSGE.get_class_longname(class)
            transforms = read(file, string(class_long) * "_revtransforms")
            transform = DSGE.parse_transform(transforms[var])

            decomp_series, transform
        end

        # Reverse transform
        transformed_decomp = DSGE.scenario_mb_reverse_transform(decomp_series, transform, :forecast)

        # Compute mean and add to DataFrame
        decomp[comp] = vec(mean(transformed_decomp, 1))
    end

    # Individual shock components are ndraws x nperiods x nshocks
    input_file = input_files[Symbol(:decompindshock, class)]
    if isfile(input_file)
        jldopen(input_file, "r") do file
            shock_indices = read(file, "shock_indices")
            for (shock, i) in shock_indices
                # Read in variable and shock: ndraws x nperiods
                decomp_series = DSGE.read_forecast_series(file, class, :decompindshock, var, shock)

                # Parse transform
                class_long = DSGE.get_class_longname(class)
                transforms = read(file, string(class_long) * "_revtransforms")
                transform = DSGE.parse_transform(transforms[var])

                # Reverse transform
                transformed_decomp = DSGE.scenario_mb_reverse_transform(decomp_series, transform, :forecast)

                # Compute mean and add to DataFrame
                decomp[shock] = vec(mean(transformed_decomp, 1))
            end
        end
    end

    return decomp
end
