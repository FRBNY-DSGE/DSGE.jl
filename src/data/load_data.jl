"""
```
load_data(m::AbstractModel; start_date = "", end_date = "", population_mnemonic = :CNP16OV
```

Loads the data for model `m` and transforms it appropriately. Returns a DataFrame.

### Inputs
- `m`: Model object

#### Keyword arguments
- `start_date`: first date to obtain data for. Default = 06/30/1959.
- `end_date`: last date to obtain data for. Default = last day of current quarter.
- `population_mnemonic`: the name of the column that will hold the
  population measure for computing per-capita values. By default, it
  is CNP16OV (Civilian Noninstitutional Population, in thousands,
  obtained via the FRED API).
"""
function load_data(m::AbstractModel; start_date="", end_date="", population_mnemonic = :CNP16OV)

    start_date = if isempty(start_date)
        Date("1959-06-30","y-m-d")
    else
        Date(start_date,"y-m-d")
    end

    end_date = if isempty(end_date)
        last_quarter_end()
    else
        Date(end_date, "y-m-d")
    end

    levels = load_levels_data(m, start_date=start_date, end_date=end_date)
    transformed = transform_data(m, levels, population_mnemonic)
    
    return transformed
end


"""
```
load_levels_data(m::AbstractModel; Date("1959-03-31","y-m-d")::Date, end_date=last_quarter_end()::Date)
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
function load_data(m::AbstractModel; start_date::Date = Date("1959-03-31","y-m-d"),
                                     end_date::Date   = last_quarter_end())


    # Load FRED data, set ois series to load
    data = load_fred_data(m; start_date=firstdayofquarter(start_date), end_date=end_date)
    set_ois_series!(m)

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
        file = inpath(m, "data", "$(string(source))_$vint.csv")
        if isfile(file)
            println("Reading data from $file...")

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
`set_ois_series!(m::AbstractModel)`

Sets the series to load from the file `ois_vint.csv` based on `n_anticipated_shocks(m)`.
"""
function set_ois_series!(m::AbstractModel)
    nant = n_anticipated_shocks(m)
    if nant > 0
        m.data_series[:ois] = [symbol("ant$i") for i in 1:nant]
    end
end
