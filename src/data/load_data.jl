"""
```
load_data(m::AbstractModel)
```

Create a DataFrame with all data series for this model, fully transformed.  

Check in `inpath(m)` for vintaged datasets corresponding to the ones in
`keys(m.data_series)`. Load the appriopriate data series (specified in
`m.data_series[source]`) for each data source. The start and end dates of the data returned
are given by the date_presample_start and date_mainsample_end model settings. Save the
resulting DataFrame to disk.

# Notes

The name of the input data file must be the same as the source string in `m.data_series`,
and those files must be located in .csv files in `inpath(m, "data")`. To accomodate growth
rates and other similar transformations, more rows of data may be downloaded than otherwise
specified by the date model settings.
"""
function load_data(m::AbstractModel)
    df = load_data_levels(m)
    df = transform_data(m, df)

    # Ensure that only appropriate rows make it into the returned DataFrame.
    start_date = get_setting(m, :date_presample_start)
    end_date   = get_setting(m, :date_mainsample_end)
    df = df[start_date .<= df[:, :date] .<= end_date, :]

    save_data(m, df)

    return df
end

function load_data_levels(m::AbstractModel)
    # Start two quarters further back than `start_date` as we need these additional
    # quarters to compute differences.
    start_date = get_setting(m, :date_presample_start) - Dates.Month(6)
    end_date = get_setting(m, :date_mainsample_end)

    # Load FRED data
    df = load_fred_data(m; start_date=firstdayofquarter(start_date), end_date=end_date)

    # Set ois series to load
    if n_anticipated_shocks(m) > 0
        m.data_series[:ois] = [symbol("ant$i") for i in 1:n_anticipated_shocks(m)]
    end

    # For each additional source, search for the file with the proper name. Open
    # it, read it in, and merge it with fred_series

    vint = data_vintage(m)

    for source in keys(m.data_series)

        # Skip FRED sources, which are handled separately
        if source == :fred
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
            df = join(df, addl_data, on=:date, kind=:outer)
        else
            # If series not found, use all NaNs
            addl_data = DataFrame(fill(NaN, (size(df,1), length(mnemonics))))
            names!(addl_data, mnemonics)
            df = hcat(df, addl_data)
            warn("$file was not found; NaNs used." )
        end
    end

    # turn NAs into NaNs
    na2nan!(df)

    sort!(df, cols = :date)

    return df
end

"""
```
save_data(m::AbstractModel, df::DataFrame)
```

Save `df` to disk as CSV. File is located in `inpath(m, "data")`. Note that this file is not
currently used for anything besides convenient visual inspection.
"""
function save_data(m::AbstractModel, df::DataFrame)
    vint = data_vintage(m)
    filename = inpath(m, "data", "data_$vint.csv")
    writetable(filename, df; nastring="NaN")
end

"""
```
df_to_matrix(m::AbstractModel, df::DataFrame)
```

Return `df`, converted to matrix of floats, and discard date column.
"""
function df_to_matrix(m::AbstractModel, df::DataFrame)
    #TODO sort columns as well according to observables
    datecol = df.colindex[:date]
    inds = setdiff(1:size(df,2), datecol)
    return convert(Matrix{Float64}, df[:, inds])
end
