function quarter_date_to_number(date::Date)
    y = Dates.year(date)
    m = Dates.month(date)
    if m == 3
        return y
    elseif m == 6
        return y + 0.25
    elseif m == 9
        return y + 0.5
    elseif m == 12
        return y + 0.75
    end
end

function quarter_number_to_date(datenum::Real)
    if datenum % 0.25 != 0
        throw(DomainError())
    end

    y = convert(Int, floor(datenum))
    q = convert(Int, (datenum % 1) / 0.25) + 1
    return quartertodate("$y-Q$q")
end

function get_date_ticks(start_date::Date, end_date::Date;
                        tick_size::Int = 5)
    dates = quarter_range(start_date, end_date)
    get_date_ticks(dates, tick_size = tick_size)
end

function get_date_ticks(dates::AbstractArray{Date, 1};
                        tick_size::Int = 5)
    datenums = map(quarter_date_to_number, dates)
    t0 = convert(Int, ceil(datenums[1] / tick_size) * tick_size)
    t1 = convert(Int, floor(datenums[end]))
    ticks = t0:tick_size:t1
    return ticks
end

function shockdec_date_to_x(date::Date, start_date::Date)
    start_x = 0.5
    quarters_diff = subtract_quarters(date, start_date)
    x = start_x + quarters_diff
    return x
end

function date_ticks!(p::Plots.Plot,
                     start_date::Date, end_date::Date,
                     tick_size::Int)
    # xlims attribute only sets finite values, e.g. (-Inf, 2) sets only the right limit
    t0 = quarter_date_to_number(start_date)
    t1 = quarter_date_to_number(end_date)

    date_ticks = get_date_ticks(start_date, end_date, tick_size = tick_size)
    xaxis!(p, xlims = (t0, t1), xtick = date_ticks)

    return nothing
end

function get_date_limit_indices(start_date::Date, end_date::Date,
                                dates::AbstractArray{Date, 1})
    start_ind = if dates[1] <= start_date <= dates[end]
        findfirst(dates, start_date)
    elseif start_date < dates[1]
        1
    else
        error("start_date $start_date cannot be after last forecast period $(dates[end])")
    end

    end_ind = if dates[1] <= end_date <= dates[end]
        findfirst(dates, end_date)
    elseif end_date > dates[end]
        length(dates)
    else
        error("end_date $end_date cannot be before first historical period $(dates[1])")
    end

    return start_ind, end_ind
end

function has_nonidentical_bands(var::Symbol, mb::MeansBands)
    nanapprox(x, y) = x ≈ y || (isnan(x) && isnan(y))

    df = mb.bands[var]
    cols = setdiff(names(df), [:date])
    for t = 1:size(df, 1)
        bandvals = convert(Matrix, df[t, cols])
        if !all(x -> nanapprox(x, mb.means[t, var]), bandvals)
            return true
        end
    end
    return false
end

function get_bands_indices(var::Symbol, history::MeansBands, forecast::MeansBands,
                           hist_inds::UnitRange{Int}, fcast_inds::UnitRange{Int})

    hist_bands  = !isempty(history)  && has_nonidentical_bands(var, history)
    fcast_bands = !isempty(forecast) && has_nonidentical_bands(var, forecast)

    if hist_bands && fcast_bands
        return hist_inds.start:fcast_inds.stop
    elseif hist_bands
        return hist_inds
    elseif fcast_bands
        return fcast_inds
    else
        return 1:0
    end
end

function plot_extension()
    be = typeof(Plots.backend())
    if be == Plots.GRBackend
        :pdf
    elseif be in [Plots.PlotlyBackend, Plots.PlotlyJSBackend]
        :html
    else
        :pdf
    end
end

function describe_series(m::AbstractModel, var::Symbol, class::Symbol;
                         detexify::Bool = false)
    res = if class in [:obs, :pseudo]
        dict = if class == :obs
            m.observable_mappings
        elseif class == :pseudo
            pseudo_measurement(m)[1]
        end
        dict[var].name
    elseif class == :states
        string(var)
    elseif class in [:shocks, :stdshocks]
        replace(string(var), r"_sh$", "")
    else
        error("Invalid class: " * string(class))
    end

    detexify ? DSGE.detexify(res) : res
end

function series_ylabel(m::AbstractModel, var::Symbol, class::Symbol;
                       untrans::Bool = false,
                       fourquarter::Bool = false)
    if untrans && fourquarter
        error("Only one of untrans or fourquarter can be true")
    end

    if class in [:obs, :pseudo]
        dict = if class == :obs
            m.observable_mappings
        elseif class == :pseudo
            pseudo_measurement(m)[1]
        end
        transform = dict[var].rev_transform

        if transform in [loggrowthtopct_annualized_percapita, loggrowthtopct_annualized]
            if untrans
                return "Q/Q Log Growth Rate"
            elseif fourquarter
                return "Percent 4Q Growth"
            else
                return "Percent Q/Q Annualized"
            end
        elseif transform in [logleveltopct_annualized_percapita, logleveltopct_annualized]
            if untrans
                return "Log Level"
            elseif fourquarter
                return "Percent 4Q Growth"
            else
                return "Percent Q/Q Annualized"
            end
        elseif transform == quartertoannual
            if untrans
                return "Percent Q/Q"
            else
                return "Percent Annualized"
            end
        elseif transform == identity
            ""
        end
    elseif class == :stdshocks
        return "Standard Deviations"
    elseif class in [:states, :shocks]
        return ""
    else
        error("Invalid class: " * string(class))
    end
end

function save_plot(p::Plots.Plot, output_file::String = "")
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkpath(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end
end