"""
```
hair_plot(var, realized, histories, forecasts; kwargs...)

hair_plot(var, realized, initial_values, forecasts; plotroot = "",
    verbose = :low, kwargs...)
```

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `realized::DataFrame`: must contain realized values of `var`
- `histories::Vector{MeansBands}` (method 1) or
  `initial_values::Vector{Float64}` (method 2): vector of either historical
  `MeansBands` or initial forecast values (i.e. s_{T|T} or y_T). Needed to
  connect the forecast hairs to the realized data line
- `forecasts::Vector{MeansBands}`

### Keyword Arguments

- `plotroot::String`: if nonempty, plots will be saved in that directory
- `verbose::Symbol`

See `?hair` for additional keyword arguments, all of which can be passed into
`hair_plot`.

### Output

- `p::Plot`
"""
function hair_plot(var::Symbol, realized::DataFrame,
                   histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                   kwargs...)

    initial_values = map(history -> history.means[end, var], histories)
    hair_plot(var, realized, initial_values, forecasts; kwargs...)
end

function hair_plot(var::Symbol, realized::DataFrame,
                   initial_values::Vector{Float64}, forecasts::Vector{MeansBands};
                   plotroot::String = "",
                   verbose::Symbol = :low,
                   kwargs...)
    # Call recipe
    p = hair(var, realized, initial_values, forecasts; kwargs...)

    # Save plot
    if !isempty(plotroot)
        output_file = joinpath(plotroot, "hairplot_" * detexify(string(var)) * "." *
                               string(plot_extension()))
        save_plot(p, output_file, verbose = verbose)
    end

    return p
end

@userplot Hair

"""
```
hair(var, realized, initial_values, forecasts;
    hist_label = \"Realized\", forecast_label = \"Forecasts\",
    hist_color = :black, forecast_color = :red, forecast_palette = :none,
    tick_size = 2)
```

User recipe called by `hair_plot`.

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `initial_values::Vector{Float64}`: vector of initial forecast values (i.e. s_{T|T} or y_T). Needed to
  connect the forecast hairs to the realized data line
- `forecasts::Vector{MeansBands}`

### Keyword Arguments

- `hist_label::String`
- `forecast_label::String`
- `forecast_color`
- `forecast_palette`: if not `:none`, the hair colors will be chosen according
  to this palette; otherwise they will all be `forecast_color`. Values
  correspond to values of the Plots attribute `color_palette` (see
  docs.juliaplots.org/latest/attributes)
- `tick_size::Int`: x-axis (time) tick size in units of years

Additionally, all Plots attributes (see docs.juliaplots.org/latest/attributes)
are supported as keyword arguments.
"""
hair

@recipe function f(hp::Hair;
                   hist_label = "Realized",
                   forecast_label = "Forecasts",
                   hist_color = :black,
                   forecast_color = :red,
                   forecast_palette = :none,
                   tick_size = 2)
    # Error checking
    if length(hp.args) != 4 || typeof(hp.args[1]) != Symbol || typeof(hp.args[2]) != DataFrame ||
        !(typeof(hp.args[3]) <: AbstractVector) ||
        !(typeof(hp.args[4]) <: AbstractVector{MeansBands})

        error("hair must be given Tuple{Symbol, DataFrame, AbstractVector, AbstractVector{MeansBands}}. Got $(typeof(hf.args))")
    end

    var, realized, initial_values, forecasts = hp.args
    if length(initial_values) != length(forecasts)
        error("Lengths of initial_values ($length(initial_values)) and forecasts ($length(forecasts)) do not match")
    end

    # Assign date ticks
    date_ticks = Base.filter(x -> Dates.month(x) == 3,            realized[!, :date])
    date_ticks = Base.filter(x -> Dates.year(x) % tick_size == 0, date_ticks)
    xticks --> (date_ticks, map(Dates.year, date_ticks))

    # Realized series
    @series begin
        seriestype := :line
        linewidth := 2
        linecolor := hist_color
        label     := hist_label

        realized[!, :date], realized[!, var]
    end

    # Forecasts
    for (initial_value, forecast) in zip(initial_values, forecasts)
        @series begin
            seriestype := :line
            linewidth  := 1
            label      := forecast == forecasts[1] ? forecast_label : ""
            if forecast_palette == :none
                linecolor := forecast_color
            else
                palette   := forecast_palette
            end

            initial_date = iterate_quarters(forecast.means[1, :date], -1)
            x = vcat(initial_date,  forecast.means[!, :date])
            y = vcat(initial_value, forecast.means[!, var])
            x, y
        end
    end
end
