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

function get_date_ticks(start_date::Date, end_date::Date;
                        tick_size::Int = 5)
    dates = DSGE.quarter_range(start_date, end_date)
    get_date_ticks(dates, tick_size = tick_size)
end

function get_date_ticks(dates::AbstractArray{Date, 1};
                        tick_size::Int = 5)
    datenums = map(quarter_date_to_number, dates)
    t0 = ceil(datenums[1] / tick_size) * tick_size
    t1 = datenums[end]
    ticks = t0:tick_size:t1
    return ticks
end

function date_ticks!(p::Plots.Plot,
                     start_date::Nullable{Date}, end_date::Nullable{Date},
                     tick_size::Int)
    # xlims attribute only sets finite values, e.g. (-Inf, 2) sets only the right limit
    t0 = isnull(start_date) ? -Inf : quarter_date_to_number(get(start_date))
    t1 = isnull(end_date)   ?  Inf : quarter_date_to_number(get(end_date))

    date_ticks = get_date_ticks(get(start_date), get(end_date), tick_size = tick_size)
    xaxis!(p, xlims = (t0, t1), xtick = date_ticks)

    return nothing
end

function get_date_limits(start_date::Nullable{Date}, end_date::Nullable{Date},
                         dates::AbstractArray{Date, 1})
    if isnull(start_date)
        start_date = Nullable(dates[1])
    end
    if isnull(end_date)
        end_date = Nullable(dates[end])
    end

    return start_date, end_date
end

function get_date_limit_indices(start_date::Nullable{Date}, end_date::Nullable{Date},
                                dates::AbstractArray{Date, 1})
    start_ind = if isnull(start_date)
        1
    else
        if dates[1] <= get(start_date) <= dates[end]
            findfirst(dates, get(start_date))
        elseif get(start_date) < dates[1]
            1
        else
            error("start_date $(get(start_date)) cannot be after last forecast period $(dates[end])")
        end
    end
    end_ind = if isnull(end_date)
        length(dates)
    else
        if dates[1] <= get(end_date) <= dates[end]
            findfirst(dates, get(end_date))
        elseif get(end_date) > dates[end]
            length(dates)
        else
            error("end_date $(get(end_date)) cannot be before first historical period $(dates[1])")
        end
    end
    return start_ind, end_ind
end

function save_plot(p::Plots.Plot, output_file::String = "")
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkdir(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end
end