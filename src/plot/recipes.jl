@userplot HistForecast

@recipe function f(hf::HistForecast;
                   start_date = hf.args[2].means[1, :date],
                   end_date = hf.args[3].means[end, :date],
                   hist_label = "History",
                   forecast_label = "Forecast",
                   hist_color = :black,
                   forecast_color = :red,
                   bands_color = :blue,
                   bands_pcts = union(which_density_bands(hf.args[2], uniquify = true),
                                      which_density_bands(hf.args[3], uniquify = true)),
                   bands_style = :fan,
                   label_bands = false,
                   transparent_bands = true,
                   tick_size = 2)
    # Error checking
    if length(hf.args) != 3 || typeof(hf.args[1]) != Symbol ||
        typeof(hf.args[2]) != MeansBands || typeof(hf.args[3]) != MeansBands

        error("histforecast must be given a Symbol and two MeansBands. Got $(typeof(hf.args))")
    end

    # Concatenate MeansBands
    var, hist, forecast = hf.args
    combined = cat(hist, forecast)
    dates = combined.means[:date]

    # Assign date ticks
    date_ticks = Base.filter(x -> start_date <= x <= end_date,    dates)
    date_ticks = Base.filter(x -> Dates.month(x) == 3,            date_ticks)
    date_ticks = Base.filter(x -> Dates.year(x) % tick_size == 0, date_ticks)
    xticks --> (date_ticks, map(Dates.year, date_ticks))

    # Bands
    sort!(bands_pcts, rev = true) # s.t. non-transparent bands will be plotted correctly
    inds = find(start_date .<= combined.bands[var][:date] .<= end_date)

    for (i, pct) in enumerate(bands_pcts)
        seriestype := :line

        x = combined.bands[var][inds, :date]
        lb = combined.bands[var][inds, Symbol(pct, " LB")]
        ub = combined.bands[var][inds, Symbol(pct, " UB")]

        if bands_style == :fan
            @series begin
                if transparent_bands
                    fillcolor := bands_color
                    fillalpha := 0.1
                else
                    if typeof(bands_color) in [Symbol, String]
                        bands_color = parse(Colorant, bands_color)
                    end
                    fillcolor := weighted_color_mean(0.1*i, bands_color, colorant"white")
                    fillalpha := 1
                end
                linealpha  := 0
                fillrange  := ub
                label      := label_bands ? "$pct Bands" : ""
                x, lb
            end
        elseif bands_style == :line
            # Lower bound
            @series begin
                linecolor := bands_color
                label     := label_bands ? "$pct LB" : ""
                x, lb
            end

            # Upper bound
            @series begin
                linecolor := bands_color
                label     := label_bands ? "$pct UB" : ""
                x, ub
            end
        else
            error("bands_style must be either :fan or :line. Got $bands_style")
        end
    end

    # Mean history
    @series begin
        seriestype :=  :line
        linewidth  --> 2
        linecolor  :=  hist_color
        label      :=  hist_label

        inds = intersect(find(start_date .<= dates .<= end_date),
                         find(hist.means[1, :date] .<= dates .<= hist.means[end, :date]))
        combined.means[inds, :date], combined.means[inds, var]
    end

    # Mean forecast
    @series begin
        seriestype :=  :line
        linewidth  --> 2
        linecolor  :=  forecast_color
        label      :=  forecast_label

        inds = intersect(find(start_date .<= dates .<= end_date),
                         find(hist.means[end, :date] .<= dates .<= forecast.means[end, :date]))
        combined.means[inds, :date], combined.means[inds, var]
    end
end