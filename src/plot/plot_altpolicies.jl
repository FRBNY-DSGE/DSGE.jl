function plot_altpolicies{T<:AbstractModel}(models::Vector{T}, vars::Vector{Symbol}, class::Symbol,
                                            cond_type::Symbol;
                                            forecast_string::String = "",
                                            altpol_string::String = "",
                                            fourquarter::Bool = false,
                                            plotroot::String = figurespath(m, "forecast"),
                                            titles::Vector{String} = String[],
                                            start_date = iterate_quarters(date_mainsample_end(m), -4),
                                            end_date = iterate_quarters(date_forecast_start(m), 20),
                                            kwargs...)
    n_altpolicies = length(models)

    # Determine output_vars
    if fourquarter
        hist_prod  = :hist4q
        fcast_prod = :forecast4q
    else
        hist_prod  = :hist
        fcast_prod = :forecast
    end

    # Read in MeansBands
    read_mb(models[1], :mode, cond_type, Symbol(hist_prod, class),
            forecast_string = forecast_string)
    mb_hists  = map(m -> read_mb(m, :mode, cond_type, Symbol(hist_prod, class),
                                 forecast_string = forecast_string),
                    models)
    mb_fcasts = map(m -> read_mb(m, :mode, cond_type, Symbol(fcast_prod, class),
                                 forecast_string = forecast_string),
                    models)

    # Determine plotting styles
    altpolicies     = map(alternative_policy, models)
    hist_labels     = fill("", n_altpolicies)
    forecast_labels = map(string, altpolicies)
    forecast_colors = Colorant[altpol.color for altpol in altpolicies]
    linestyles      = map(altpol -> altpol.linestyle, altpolicies)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(models[1], var, class, detexify = detexify_title), vars)
    end

    # Temporarily set models[1] altpolicy key to altpol_string for get_forecast_filename
    if n_altpolicies != 1
        if isempty(altpol_string)
            error("Must provide nonempty altpol_key if plotting multiple alternative policies")
        else
            altpol = get_setting(models[1], :alternative_policy)
            old_key = altpol.key
            altpol.key = Symbol(altpol_string)
        end
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        output_file = if isempty(plotroot)
            ""
        else
            get_forecast_filename(plotroot, filestring_base(models[1]), :mode, cond_type,
                                                            Symbol("altpol", fcast_prod, "_", detexify(var)),
                                  forecast_string = forecast_string,
                                  fileformat = plot_extension())
        end

        p = plot_history_and_forecast(var, mb_hists, mb_fcasts;
                                      start_date = start_date,
                                      end_date = end_date,
                                      output_file = output_file,
                                      title = title,
                                      hist_label = hist_labels,
                                      forecast_label = forecast_labels,
                                      forecast_color = forecast_colors,
                                      linestyle = linestyles,
                                      tick_size = 2,
                                      kwargs...)
    end

    # Reset models[1] altpolicy key to original value
    if n_altpolicies != 1
        altpol.key = old_key
    end

    return plots
end