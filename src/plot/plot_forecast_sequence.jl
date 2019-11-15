"""
```
plot_forecast_sequence(models, input_types, cond_types,
                       model_realized, input_type_realized, cond_type_realized,
                       class, forecast_var;
                       forecast_strings = fill("", length(models)),
                       forecast_string_realized = "",
                       title = string(forecast_var),
                       start_date::Date = quartertodate("2014-Q4"),
                       end_date::Date = quartertodate("2021-Q4"),
                       forecast_display_length = 4,
                       filepath::String = "")
```
Plot a sequence of forecasts (pulled from the tuple (model, input_type, cond_type, class, forecast_string)) on the same figure along with a realized series given by (model_realized, ...).
This is useful for seeing the evolution of forecasts with respect to the true realization of a series over time, e.g. is the DSGE model constantly under-forecasting realized inflation?

### Arguments
- `models::AbstractVector`: A vector of model objects with the correct vintage/saveroot settings from which to load the sequence of forecast MeansBands.
- `input_types::Vector{Symbol}`: A vector of the input types of each of the forecasts that were run.
- `cond_types::Vector{Symbol}`: A vector of the cond types of each of the forecasts that were run.
- `model_realized::AbstractDSGEModel`: The model object from which to load the historical realization/smoothed history of the series to compare to the forecasts.
- `input_type_realized::Symbol`: The input type from which the historical realization/smoothed history was generated.
- `cond_type_realized::Symbol`: The cond type from which the historical realization/smoothed history was generated.
- `class::Symbol`: The class, either :obs or :pseudo, that `forecast_var` belongs to.
- `forecast_var::Symbol`: The forecasted variable that is being plotted, e.g. :obs_gdp, or :NaturalRate.

### Keyword Arguments
- `forecast_strings::Vector{String}`: A vector of the forecast_strings of each of the forecasts that were run.
- `forecast_string_realized::String`: The forecast_string from which the historical realization/smoothed history was generated.
- `title::String`: The title.
- `start_date::Date`: The start date.
- `end_date::Date`: The end date.
- `forecast_display_length::String`: The number of periods of forecast that you want plotted.
- `filepath::String`: The filepath that you want to save the figure to.
"""
function plot_forecast_sequence(models::AbstractVector,
                                input_types::Vector{Symbol},
                                cond_types::Vector{Symbol},
                                model_realized::AbstractDSGEModel,
                                input_type_realized::Symbol,
                                cond_type_realized::Symbol,
                                class::Symbol,
                                forecast_var::Symbol;
                                forecast_strings::Vector{String} = fill("", length(models)),
                                forecast_string_realized::String = "",
                                title::String = string(forecast_var),
                                start_date::Date = quartertodate("2014-Q4"),
                                end_date::Date = quartertodate("2021-Q4"),
                                forecast_display_length::Int = 4,
                                plot_to_realized_end_date::Bool = false,
                                realized_linecolor::Symbol = :gray,
                                realized_linestyle::Symbol = :dash,
                                constant_to_plot::Float64 = Inf,
                                constant_linecolor::Symbol = :black,
                                ylims::Tuple = Tuple{}(),
                                filepath::String = "", kwargs...)
    n_models = length(models)
    @assert length(input_types) == length(cond_types) == length(forecast_strings) == n_models

    mb_realized = read_mb(model_realized, input_type_realized,
                          cond_type_realized, Symbol("hist", class),
                          forecast_string = forecast_string_realized)

    # Read in Forecast MBs
    mbs_forecast     = Vector{MeansBands}(undef, n_models)
    for i in 1:n_models
        mbs_forecast[i] = read_mb(models[i], input_types[i], cond_types[i],
                                  Symbol("forecast", class),
                                  forecast_string = forecast_strings[i])
    end

    plot_forecast_sequence(mbs_forecast, mb_realized, forecast_var;
                           title = title, start_date = start_date,
                           forecast_display_length = forecast_display_length,
                           plot_to_realized_end_date = plot_to_realized_end_date,
                           realized_linecolor = realized_linecolor, realized_linestyle = realized_linestyle,
                           constant_to_plot = constant_to_plot, constant_linecolor = constant_linecolor,
                           ylims = ylims, filepath = filepath)
end

function plot_forecast_sequence(mbs_forecast::Vector{MeansBands},
                                mb_realized::MeansBands,
                                forecast_var::Symbol;
                                title::String = string(forecast_var),
                                start_date::Date = quartertodate("2014-Q4"),
                                end_date::Date = quartertodate("2021-Q4"),
                                forecast_display_length::Int = 4,
                                plot_to_realized_end_date::Bool = false,
                                realized_linecolor::Symbol = :gray,
                                realized_linestyle::Symbol = :dash,
                                constant_to_plot::Float64 = Inf,
                                constant_linecolor::Symbol = :black,
                                ylims::Tuple = Tuple{}(),
                                filepath::String = "", kwargs...)
    n_forecasts = length(mbs_forecast)

    df_realized = mb_realized.means

    dfs_forecast = Vector{DataFrame}(undef, n_forecasts)
    for i in 1:n_forecasts
        dfs_forecast[i] = mbs_forecast[i].means
    end

    plot_forecast_sequence(dfs_forecast, df_realized, forecast_var;
                           title = title, start_date = start_date,
                           end_date = end_date, forecast_display_length = forecast_display_length,
                           plot_to_realized_end_date = plot_to_realized_end_date,
                           realized_linecolor = realized_linecolor, realized_linestyle = realized_linestyle,
                           constant_to_plot = constant_to_plot, constant_linecolor = constant_linecolor,
                           ylims = ylims, filepath = filepath)
end

function plot_forecast_sequence(dfs_forecast::Vector{DataFrame},
                                df_realized::DataFrame,
                                forecast_var::Symbol;
                                title::String = string(forecast_var),
                                start_date::Date = quartertodate("2014-Q4"),
                                end_date::Date = quartertodate("2021-Q4"),
                                forecast_display_length::Int = 4,
                                plot_to_realized_end_date::Bool = false,
                                realized_linecolor::Symbol = :gray,
                                realized_linestyle::Symbol = :dash,
                                constant_to_plot::Float64 = Inf,
                                constant_linecolor::Symbol = :black,
                                ylims::Tuple = Tuple{}(),
                                filepath::String = "")
    n_forecasts = length(dfs_forecast)

    # Construct realized series date indices
    realized_date_range_raw = Base.filter(x -> start_date <= x <= end_date, df_realized[!,:date])
    realized_date_range = map(quarter_date_to_number, realized_date_range_raw)
    date_inds  = findall(start_date .<= df_realized[!,:date] .<= end_date)

    # Plot realized
    p = plot(realized_date_range, df_realized[!,forecast_var][date_inds],
             label = "Realized", title = title, linecolor = realized_linecolor,
             linestyle = realized_linestyle, ylims = ylims)

    if !isinf(constant_to_plot)
        constant_start_ind   = realized_date_range[1]
        constant_end_ind     = realized_date_range[end]
    end

    # Plot forecasts
    for (i, df) in zip(1:n_forecasts, dfs_forecast)

        if plot_to_realized_end_date
            realized_end_date = df_realized[:date][end]

            # If `df[:date][1]` (the first forecast quarter)
            # is less than `forecast_display_length` away from the
            # realized_end_date, then plot the forecast up to
            # `forecast_display_length` number of periods.
            # Else, plot the forecast up to realized_end_date.
            # E.g. For forecasts that originate close to (or past) realized_end_date,
            # we still want to see a full forecast, but for forecasts that originate
            # far earlier than the current realized_end_date, we want to see their
            # forecast up to the realized_end_date.
            if subtract_quarters(realized_end_date, df[:date][1]) < forecast_display_length
                # Construct forecast date indices
                forecast_date_range_raw = Base.filter(x -> start_date <= x <= end_date,
                                                      dfs_forecast[i][:date])
                forecast_date_range = map(quarter_date_to_number, forecast_date_range_raw)

                # Truncate
                forecast_date_range = forecast_date_range[1:min(forecast_display_length, length(forecast_date_range))]
                date_inds = findall(start_date .<= dfs_forecast[i][:date] .<= end_date)
                date_inds = date_inds[1:min(forecast_display_length, length(forecast_date_range))]
            else
                # Construct forecast date indices
                forecast_date_range_raw = Base.filter(x -> start_date <= x <= realized_end_date,
                                                      dfs_forecast[i][:date])
                forecast_date_range = map(quarter_date_to_number, forecast_date_range_raw)

                date_inds = findall(start_date .<= dfs_forecast[i][:date] .<= realized_end_date)
            end
        else
            # Construct forecast date indices
            forecast_date_range_raw = Base.filter(x -> start_date <= x <= end_date,
                                                  dfs_forecast[i][!,:date])
            forecast_date_range = map(quarter_date_to_number, forecast_date_range_raw)

            # Truncate
            forecast_date_range = forecast_date_range[1:min(forecast_display_length, length(forecast_date_range))]
            date_inds = date_inds[1:min(forecast_display_length, length(date_inds))]
        end

        if i == 1
            label = "Forecast"
        else
            label = ""
        end

        plot!(p, forecast_date_range, dfs_forecast[i][!,forecast_var][date_inds],
              label = label, linecolor = :red, legend = :bottomright)

        if !isinf(constant_to_plot)
            constant_end_ind = max(constant_end_ind, forecast_date_range[end])
        end
    end

    # Plot the constant
    if !isinf(constant_to_plot)
        plot!(p, constant_start_ind:constant_end_ind,
              fill(constant_to_plot, length(constant_start_ind:constant_end_ind)),
              linecolor = constant_linecolor, label = "")
    end

    if !isempty(filepath)
        savefig(p, filepath)
        println("Saved $filepath")
    end

    return p
end
