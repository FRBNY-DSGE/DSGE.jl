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
    return DSGE.quartertodate("$y-Q$q")
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
    quarters_diff = DSGE.subtract_quarters(date, start_date)
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

function get_date_limits(nullable_start::Nullable{Date}, nullable_end::Nullable{Date},
                         dates::AbstractArray{Date, 1})

    start_date = isnull(nullable_start) ? dates[1] : get(nullable_start)
    end_date = isnull(nullable_end) ? dates[end] : get(nullable_end)

    return start_date, end_date
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
    df = mb.bands[var]
    cols = setdiff(names(df), [:date])
    for t = 1:size(df, 1)
        bandvals = convert(Matrix, df[t, cols])
        if !all(x -> x ≈ mb.means[t, var], bandvals)
            return true
        end
    end
    return false
end

function get_bands_indices(var::Symbol, history::MeansBands, forecast::MeansBands,
                           hist_inds::UnitRange{Int}, fcast_inds::UnitRange{Int})

    hist_bands  = has_nonidentical_bands(var, history)
    fcast_bands = has_nonidentical_bands(var, forecast)

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

function save_plot(p::Plots.Plot, output_file::String = "")
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkdir(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end
end