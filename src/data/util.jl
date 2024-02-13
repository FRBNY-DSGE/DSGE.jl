"""
`prev_quarter(q::TimeType = now())`

Returns Date identifying last day of the previous quarter
"""
function prev_quarter(q::TimeType = now())
    Dates.Date(lastdayofquarter(firstdayofquarter(q)-Dates.Day(1)))
end

"""
`next_quarter(q::TimeType = now())`

Returns Date identifying last day of the next quarter
"""
function next_quarter(q::TimeType = now())
    Dates.Date(lastdayofquarter(lastdayofquarter(q)+Dates.Day(1)))
end

"""
`get_quarter_ends(start_date::Date,end_date::Date)`

Returns an Array of quarter end dates between `start_date` and `end_date`.
"""
function get_quarter_ends(start_date::Dates.Date,end_date::Dates.Date)
    map(lastdayofquarter, collect(start_date:Dates.Month(3):end_date))
end

"""
`quartertodate(string::String)`

Convert `string` in the form "YYqX", "YYYYqX", or "YYYY-qX" to a Date of the end of
the indicated quarter. "X" is in `{1,2,3,4}` and the case of "q" is ignored.
"""
function quartertodate(string::String)
    if occursin(r"^[0-9]{2}[qQ][1-4]$", string)
        year = "20"*string[1:2]
        quarter = string[end]
    elseif occursin(r"^[0-9]{4}[qQ][1-4]$", string)
        year = string[1:4]
        quarter = string[end]
    elseif occursin(r"^[0-9]{4}-[qQ][1-4]$", string)
        year = string[1:4]
        quarter = string[end]
    else
        throw(Meta.ParseError("Invalid format: $string"))
    end

    year = parse(Int, year)
    quarter = parse(Int, quarter)
    month = 3*quarter
    day = 1

    return lastdayofquarter(Date(year, month, day))
end

"""
`datetoquarter(date::Date)`

Convert `string` in the form "YYqX", "YYYYqX", or "YYYY-qX" to a Date of the end of
the indicated quarter. "X" is in `{1,2,3,4}` and the case of "q" is ignored.

Return an integer from the set `{1,2,3,4}`, corresponding to one of the quarters in a year given a Date object.
"""
function datetoquarter(date::Dates.Date)
    month = Dates.month(date)
    if month in 1:3
        return 1
    elseif month in 4:6
        return 2
    elseif month in 7:9
        return 3
    elseif month in 10:12
        return 4
    else
        throw("Must provide a date object with a valid month (in 1:12)")
    end
end

"""
```
function vinttodate(vint)
```

Return the string given by data_vintage(m), which is in the format YYYYMMDD, to a Date object.
"""
function vinttodate(vint::String)
    year  = Meta.parse("20"*vint[1:2])
    month = Meta.parse(vint[3:4])
    day   = Meta.parse(vint[5:6])
    return Date(year, month, day)
end

"""
`subtract_quarters(t1::Date, t0::Date)`

Compute the number of quarters between t1 and t0, including t0 and excluding t1.
"""
function subtract_quarters(t1::Dates.Date, t0::Dates.Date)
    days = t1 - t0
    quarters = round(days.value / 365.25 * 4.0)
    return convert(Int, quarters)
end


"""
`format_dates!(col, df)`

Change column `col` of dates in `df` from String to Date, and map any dates given in the
interior of a quarter to the last day of the quarter.
"""
function format_dates!(col::Symbol, df::DataFrame)
    df[!,col] = Dates.Date.(df[!,col])
    map!(lastdayofquarter, df[!,col], df[!,col])
end

"""
```
missing2nan(a::Array)
```

Convert all elements of Union{X, Missing.Missing} or Missing.Missing to type Float64.
"""
function missing2nan(a::Array)
    a_new = Array{Float64}(undef, size(a))
    for i in eachindex(a)
        a_new[i] = ismissing(a[i]) ? NaN : a[i]
    end
    return a_new
end

function missing2nan!(df::DataFrame)
    for col in propertynames(df)
        if col != :date
            df[!, col] = map(x -> ismissing(x) ? NaN : x, df[!, col])
            df[!, col] = convert(Vector{Float64}, df[!, col])
        end
    end
end

# Temp function to ensure proper transition into 0.7
function nan2missing!(df::DataFrame)
    for col in propertynames(df)
        if col != :date
            df[!,col] = convert(Vector{Union{Missing, Float64}}, df[!,col])
            df[!,col] = replace(x -> isnan(x) ? missing : x, df[!,col])
        end
    end
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(df::DataFrame)
    for col in propertynames(df)
        if typeof(df[!,col])==Vector{Date}
            nothing
        else
            df[ismissing.(df[!,col]), col] .= NaN
        end
    end
end

"""
```
na2nan!(df::Array)
```

Convert all NAs in an Array to NaNs.
"""
function na2nan!(v::Array)
    for i = 1:length(v)
        v[i] = isnan(v[i]) ?  NaN : v[i]
    end
end

"""
```
missing_cond_vars!(m, df; cond_type = :none, check_empty_columns = true)
```

Make conditional period variables not in `cond_semi_names(m)` or
`cond_full_names(m)` missing if necessary.
"""
function missing_cond_vars!(m::AbstractDSGEModel, df::DataFrame; cond_type::Symbol = :none,
                            check_empty_columns::Bool = true)
    if cond_type in [:semi, :full]
        # Get appropriate
        cond_names = if cond_type == :semi
            cond_semi_names(m)
        elseif cond_type == :full
            cond_full_names(m)
        end

        # Make non-conditional variables missing
        cond_names_missing = setdiff(propertynames(df), [cond_names; :date])
        for var_name in cond_names_missing
            df[!,var_name] = convert(Vector{Union{Missing, eltype(df[!,var_name])}}, df[!,var_name])
            df[df[!,:date] .>= date_forecast_start(m), var_name] .= missing
        end

        # Throw an error if any conditional variables are missing
        missing_vars = Vector{Symbol}(undef,0)
        for var in cond_names
            if any(ismissing.(df[df[!,:date] .>= date_forecast_start(m), var]))
                if check_empty_columns
                    push!(missing_vars, var)
                else
                    @warn "Missing some conditional observations for " * string(var)
                end
            end
        end
        if check_empty_columns && !isempty(missing_vars)
            # parse missing_vars
            to_print = Vector{String}(undef,0)
            for var in missing_vars
                push!(to_print, string(var) * ", ")
            end
            as_str = join(to_print)
            error("Column(s) $(as_str[1:end-2]) are missing conditional observations.")
        end
    end
end

"""
```
get_data_filename(m, cond_type)
```

Returns the data file for `m`, which depends on `data_vintage(m)`, and if
`cond_type in [:semi, :full]`, also on `cond_vintage(m)` and `cond_id(m)`.
"""
function get_data_filename(m::AbstractDSGEModel, cond_type::Symbol)
    filestrings = ["data"]

    # If writing conditional data, append conditional vintage and ID to filename
    if cond_type in [:semi, :full]
        push!(filestrings, "cdid=" * lpad(string(cond_id(m)), 2, string(0)))
        push!(filestrings, "cdvt=" * cond_vintage(m))
    end

    push!(filestrings, "dsid=" * lpad(string(data_id(m)), 2, string(0)))
    push!(filestrings, "vint=" * data_vintage(m))
    filename = join(filestrings, "_")

    return inpath(m, "data", filename * ".csv")
end

"""
```
iterate_quarters(start::Date, quarters::Int)
```

Returns the date corresponding to `start` + `quarters` quarters.

### Inputs
- `start`: starting date
- `quarters`: number of quarters to iterate forward or backward
"""
function iterate_quarters(start::Dates.Date, quarters::Int)

    next = start
    if quarters < 0
        for n = 1:-quarters
            next = Dates.toprev(next) do x
                Dates.lastdayofquarter(x) == x
            end
        end
    elseif quarters > 0
        for n = 1:quarters
            next = Dates.tonext(next) do x
                Dates.lastdayofquarter(x) == x
            end
        end
    end

    next
end

"""
```
reconcile_column_names(a::DataFrame, b::DataFrame)
```

adds columns of missings to a and b so that both have
the same set of column names.
"""
function reconcile_column_names(a::DataFrame, b::DataFrame)
    new_a_cols = setdiff(propertynames(b), propertynames(a))
    new_b_cols = setdiff(propertynames(a), propertynames(b))
    for col in new_a_cols
        a[!, col] = fill(missing, size(a, 1))
    end
    for col in new_b_cols
        b[!, col] = fill(missing, size(b, 1))
    end
    return a, b
end

"""
```
datetoymdvec(dt)
```
converts a Date to a vector/matrix holding the year, month, and date.
"""
function datetoymdvec(dt::Date)
    y = Year(dt).value
    m = Month(dt).value
    d = Day(dt).value
    return [y, m, d]
end

function datetoymdvec(dt::Vector{Date})
    return permutedims(hcat(map(x -> datetoymdvec(x), dt)...), [2, 1])
end

"""
```
quartertofloats(dt)
```
converts a Date to a floating point number based on the quarter
"""
function quartertofloats(dt::Date)
    dtvec = datetoymdvec(dt)
    out = dtvec[1]
    if dtvec[2] == 3
        return out
    elseif dtvec[2] == 6
        return out + .25
    elseif dtvec[3] == 9
        return out + .5
    else
        return out + .75
    end
end

function quartertofloats(dt::Vector{Date})
    return map(x -> quartertofloats(x), dt)
end
