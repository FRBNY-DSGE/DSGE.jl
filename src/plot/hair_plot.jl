"""
```
hair_plot(var, df, histories, forecasts; kwargs...)

hair_plot(var, df, initial_values, forecasts; output_file = "", hist_label = \"Realized\",
    forecast_label = \"Forecasts\", forecast_palette = Symbol(), forecast_color = :red,
    legend = :best)
```

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `df::DataFrame`: must contain realized values of `var`
- `histories::Vector{MeansBands}` (method 1) or
  `initial_values::Vector{Float64}` (method 2): vector of either historical
  `MeansBands` or initial forecast values (i.e. s_{T|T} or y_T). Needed to
  connect the forecast hairs to the realized data line
- `forecasts::Vector{MeansBands}`

### Keyword Arguments

- `output_file::String`: if specified, plot will be saved there as a PDF
- `hist_label::String`
- `forecast_label::String`
- `forecast_palette::Symbol`: if specified, the hair colors will be chosen
  according to this palette. Otherwise they will all be `forecast_color`
- `forecast_color::Colorant`
- `legend`

### Output

- `p::Plot`
"""
function hair_plot(var::Symbol, df::DataFrame,
                   histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                   kwargs...)

    initial_values = map(history -> history.means[end, var], histories)
    hair_plot(var, df, initial_values, forecasts; kwargs...)
end

function hair_plot(var::Symbol, df::DataFrame,
                   initial_values::Vector{Float64}, forecasts::Vector{MeansBands};
                   output_file::String = "",
                   hist_label::String = "Realized",
                   forecast_label::String = "Forecasts",
                   forecast_palette::Symbol = Symbol(),
                   forecast_color::Colorant = colorant"red",
                   legend = :best)
    # Dates
    datenums   = map(quarter_date_to_number, df[:date])
    date_ticks = get_date_ticks(df[:date])

    # Initialize plot
    p = Plots.plot(xtick = date_ticks, legend = legend)

    # Plot realized (transformed) series
    plot!(p, datenums, df[var], label = hist_label, linewidth = 2, linecolor = :black)

    # Plot each forecast
    for (initial_value, forecast) in zip(initial_values, forecasts)
        date_0 = DSGE.iterate_quarters(forecast.means[1, :date], -1)
        dates = vcat([date_0], forecast.means[:date])
        datenums = map(quarter_date_to_number, dates)

        series = vcat([initial_value], forecast.means[var])

        label = forecast == forecasts[1] ? forecast_label : ""
        if forecast_palette == Symbol()
            plot!(p, datenums, series, label = label, linewidth = 1, linecolor = forecast_color)
        else
            plot!(p, datenums, series, label = label, linewidth = 1, palette = forecast_palette)
        end
    end

    # Save if `output_file` provided
    save_plot(p, output_file)

    return p
end