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

function get_date_ticks(dates::AbstractArray{Date, 1})
    datenums = map(quarter_date_to_number, dates)
    t0 = ceil(datenums[start_ind] / 5) * 5
    t1 = datenums[end]
    ticks = t0:5:t1
    return ticks
end