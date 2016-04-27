"""
```
load_data(m::AbstractModel)
```

Checks in `inpath(m)` for vintaged datasets corresponding to the ones in
`keys(m.data_series)`. Loads the appriopriate data series (specified in
`m.data_series[source]`) for each data source, and returns a consolidated DataFrame (merged
on date). The start and end date of the data downloaded are given by the
date_presample_start and date_mainsample_end model settings.  

# Notes
The name of the input data file must be the same as the source string in `m.data_series`,
and those files must be located in .csv files in `inpath(m, "data")`.
"""
function load_data(m::AbstractModel)

    # Start two quarters further back than `start_date` as we need these additional
    # quarters to compute differences.
    start_date = get_setting(m, :date_presample_start) - Dates.Month(6)
    end_date = get_setting(m, :date_mainsample_end)

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
            # If series not found, use all NaNs
            addl_data = DataFrame(fill(NaN, (size(data,1), length(mnemonics))))
            names!(addl_data, mnemonics)
            data = hcat(data, addl_data)
            warn("$file was not found; NaNs used." )
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
