"""
```
load_data(m::AbstractModel; Date("1959-03-31","y-m-d")::Date, end_date=last_quarter_end()::Date)
```

Checks in `inpath(m)` for vintaged datasets corresponding to the ones in
`keys(m.data_series)`. Loads the appriopriate data series (specified in
`m.data_series[source]`) for each data source, and returns a consolidated DataFrame (merged
on date).

# Arguments
- `m::AbstractModel`: model object
- `start_date`: start date of entire sample
- `end_date`: end date of entire sample

# Notes
The name of the input data file must be the same as the source string in `m.data_series`,
and those files must be located in .csv files in `inpath(m, "data")`.
"""
function load_data(m::AbstractModel; start_date=Date("1959-03-31","y-m-d")::Date,
                                     end_date=last_quarter_end()::Date)

    
    # Load FRED data, set ois series to load, and set which longrate series to load
    data = load_fred_data(m; start_date=firstdayofquarter(start_date), end_date=end_date)
    set_ois_series!(m)
    set_longrate_series!(m)
    
    # For each additional source, search for the file with the proper name. Open
    # it, read it in, and merge it with fred_series

    vint = data_vintage(m)

    for source in keys(m.data_series)

        # Skip FRED sources, which are handled separately
        if isequal(source, :fred)
            continue
        end

        # Check that this source is actually used
        mnemonics = m.data_series[source]
        if isempty(mnemonics)
            warn("No series were specified from $(string(source))")
            continue
        end

        # Read and merge data from this source
        file = inpath(m, "data", "$(string(source))_$vint.txt")
        if isfile(file)

            # Read in dataset and check that the file contains data for the proper dates
            addl_data = readtable(file)

            # Convert dates from strings to quarter-end dates for date arithmetic
            format_dates!(:date, addl_data)

            # Warn on sources with incomplete data; missing data will be replaced with NaN
            # during merge.
            if !in(lastdayofquarter(start_date), addl_data[:date]) ||
               !in(lastdayofquarter(end_date), addl_data[:date])
                warn("$file does not contain the entire date range specified...you may want to update your data.")
            end

            # Make sure each mnemonic that was specified is present
            for series in mnemonics
                if !in(series, names(addl_data))
                    error("$(string(series)) is missing from $file.")
                end
            end

            # Extract just the columns and rows of the dataset we want, and merge them with
            # data
            cols = [:date; mnemonics]
            rows = start_date .<= addl_data[:date] .<= end_date

            addl_data = addl_data[rows, cols]
            data = join(data, addl_data, on=:date, kind=:outer)
        else
            error("$file was not found." )
        end
    end

    # turn NAs into NaNs
    na2nan!(data)

    return sort!(data, cols = :date)
end

"""
`set_longrate_series!(m::AbstractModel)`

    Sets the series to load from the file `longrate_vint.txt` based on whether we want to subtract term premia or not.
"""
function set_longrate_series!(m::AbstractModel)
    m.data_series[:longrate], m.data_transforms[:longrate] = if adjust_longrate(m)

        fn = function (levels)
            annualtoquarter(levels[:FYCZZA] - levels[:FTPZAC])
        end

        ["FYCZZA", "FTPZAC"], fn
    else
        fn = function (levels)
            annualtoquarter(levels[:FYCZZA])
        end
        
        ["FYCZZA"], fn
    end
end

"""
`set_ois_series!(m::AbstractModel)`

Sets the series to load from the file `ois_vint.txt` based on `n_anticipated_shocks(m)`.
"""
function set_ois_series!(m::AbstractModel)
    nant = n_anticipated_shocks(m)
    if nant > 0
        m.data_series[:ois] = [symbol("ant$i") for i in 1:nant]
    end
end


"""
`last_quarter_end()`

Returns the first day of the previous quarter
"""
function last_quarter_end()
    cqs = Dates.firstdayofquarter(now())     # current quarter start
    lqe = cqs - Dates.Day(1)                 # last quarter end

    Date(lqe)
end

"""
Returns a DataArray of quarter start dates between `start_date` and `end_date`.
"""
function get_quarter_ends(start_date,end_date)
    dr = start_date:end_date

    dates = recur(dr) do x
        lastdayofquarter(x) == x
    end
end

"""
Converts a collection of strings in "y-m-d" format to Dates.
"""
function stringstodates(stringarray)
    n = length(stringarray)
    dates = Vector{Date}(n)
    for i = 1:n
        dates[i] = Date(stringarray[i], "y-m-d")
    end

    dates
end

"""
`format_dates!(col, df)`

Change column `col` of dates in `df` from String to Date, and map any dates given in the
interior of a quarter to the last day of the quarter.
"""
function format_dates!(col, df)
    df[col] = stringstodates(df[col])
    map!(x->lastdayofquarter(x), df[col], df[col])
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(df::DataFrame)
    for col in names(df)
        df[isna(df[col]), col] = NaN
    end
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(v::DataArray)
    for i = 1:length(v)
        v[i] = isna(v[i]) ?  NaN : v[i]
    end
end
