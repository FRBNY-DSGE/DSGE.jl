"""
```
plot_shock_decomposition(m, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_shock_decomposition(m, vars, class, input_type, cond_type;
    forecast_string = "", groups = shock_groupings(m),
    plotroot::String = figurespath(m, \"forecast\"), titles = [],
    kwargs...)
```

Plot shock decomposition(s) for `var` or `vars`.

### Inputs

- `m::AbstractModel`
- `var::Symbol` or `vars::Vector{Symbol}`: variable(s) whose shock decomposition
  is to be plotted, e.g. `:obs_gdp` or `[:obs_gdp, :obs_nominalrate]`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

### Keyword Arguments

- `forecast_string::String`
- `groups::Vector{ShockGroup}`
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`
- `fileformat::Symbol`

See `?shockdec` for additional keyword arguments, all of which can be passed
into `plot_history_and_forecast`.

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_shock_decomposition(m::AbstractModel, var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  title = "", fileformat = plot_extension(),
                                  kwargs...)

    plots = plot_shock_decomposition(m, [var], class, input_type, cond_type;
                                     titles = isempty(title) ? String[] : [title],
                                     fileformat = fileformat, kwargs...)
    return plots[var]
end

function plot_shock_decomposition(m::AbstractModel, vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  groups::Vector{ShockGroup} = shock_groupings(m),
                                  plotroot::String = figurespath(m, "forecast"),
                                  titles::Vector{String} = String[],
                                  verbose::Symbol = :low,
                                  fileformat = plot_extension(),
                                  kwargs...)
    # Read in MeansBands
    output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
    mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string),
              output_vars)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(m, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Call recipe
        plots[var] = shockdec(var, mbs..., groups;
                              ylabel = series_ylabel(m, var, class) * "\n(deviations from mean)",
                              title = title, kwargs...)

        # Save plot
        if !isempty(plotroot)
            output_file = get_forecast_filename(plotroot, filestring_base(m), input_type, cond_type,
                                                Symbol("shockdec_", detexify(var)),
                                                forecast_string = forecast_string,
                                                fileformat = fileformat)
            save_plot(plots[var], output_file, verbose = verbose)
        end
    end
    return plots
end

@userplot Shockdec

"""
```
shockdec(var, shockdec, trend, dettrend, hist, forecast, groups;
    start_date = shockdec.means[1, :date],
    end_date = shockdec.means[end, :date],
    hist_label = \"Detrended History\",
    forecast_label = \"Detrended Forecast\",
    hist_color = :black, forecast_color = :red, tick_size = 5)
```

User recipe called by `plot_shock_decomposition`.

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `shockdec::MeansBands`
- `trend::MeansBands`
- `dettrend::MeansBands`
- `hist::MeansBands`
- `forecast::MeansBands`
- `groups::Vector{ShockGroup}`

### Keyword Arguments

- `start_date::Date`
- `end_date::Date`
- `hist_label::String`
- `forecast_label::String`
- `hist_color`
- `forecast_color`
- `tick_size::Int`: x-axis (time) tick size in units of years

Additionally, all Plots attributes (see docs.juliaplots.org/latest/attributes)
are supported as keyword arguments.
"""
shockdec

@recipe function f(sd::Shockdec;
                   start_date = sd.args[2].means[1, :date],
                   end_date = sd.args[2].means[end, :date],
                   hist_label = "Detrended History",
                   forecast_label = "Detrended Forecast",
                   hist_color = :black,
                   forecast_color = :red,
                   tick_size = 5)
    # Error checking
    if length(sd.args) != 7 || typeof(sd.args[1]) != Symbol ||
        typeof(sd.args[2]) != MeansBands || typeof(sd.args[3]) != MeansBands ||
        typeof(sd.args[4]) != MeansBands || typeof(sd.args[5]) != MeansBands ||
        typeof(sd.args[6]) != MeansBands || typeof(sd.args[7]) != Vector{ShockGroup}

        error("shockdec must be given a Symbol, five MeansBands, and a Vector{ShockGroup}. Got $(typeof(sd.args))")
    end

    var, shockdec, trend, dettrend, hist, forecast, groups = sd.args

    # Construct DataFrame with detrended mean, deterministic trend, and all shocks
    df = prepare_means_table_shockdec(shockdec, trend, dettrend, var,
                                      mb_hist = hist, mb_forecast = forecast,
                                      detexify_shocks = false,
                                      groups = groups)
    dates = df[:date]
    xnums = (1:length(dates)) - 0.5

    # Assign date ticks
    inds = intersect(find(x -> start_date .<= x .<= end_date,  dates),
                     find(x -> Dates.month(x) == 3,            dates),
                     find(x -> Dates.year(x) % tick_size == 0, dates))
    xticks --> (xnums[inds], map(Dates.year, dates[inds]))

    # Set date axis limits
    x0 = xnums[findfirst(dates .== start_date)]
    x1 = xnums[findfirst(dates .== end_date)]
    xlims := (x0, x1)

    # Shock contributions
    @series begin
        labels    = map(x -> x.name,  groups)
        cat_names = map(Symbol, labels)
        colors    = map(x -> x.color, groups)
        ngroups   = length(groups)

        labels     --> reshape(labels, 1, ngroups)
        color      --> reshape(colors, 1, ngroups)
        linealpha  --> 0
        bar_width  --> 1
        legendfont --> Plots.Font("sans-serif", 5, :hcenter, :vcenter, 0.0, colorant"black")

        x = df[:date]
        y = convert(Matrix{Float64}, df[cat_names])
        StatPlots.GroupedBar((x, y))
    end

    seriestype := :line
    linewidth  := 2

    # Detrended mean history
    @series begin
        linecolor := hist_color
        label     := hist_label

        inds = find(hist.means[1, :date] .<= dates .<= hist.means[end, :date])
        x = xnums[inds]
        y = convert(Vector{Float64}, df[inds, :detrendedMean])
        x, y
    end

    # Detrended mean forecast
    @series begin
        linecolor := forecast_color
        label     := forecast_label

        inds = find(hist.means[end, :date] .<= dates .<= forecast.means[end, :date])
        x = xnums[inds]
        y = convert(Vector{Float64}, df[inds, :detrendedMean])
        x, y
    end
end