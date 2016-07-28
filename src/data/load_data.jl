"""
```
load_data(m::AbstractModel; try_disk::Bool = true, verbose::Symbol = :low)
```

Create a DataFrame with all data series for this model, fully transformed.  

First, check the disk to see if a valid dataset is already stored in `inpath(m, \"data\")`. A
dataset is valid if every series in `m.data_transforms` is present and the entire sample is
contained (from `date_presample_start` to `date_mainsample_end`. If no valid dataset is
already stored, the dataset will be recreated. This check can be eliminated by passing
`try_disk=false`.

If the dataset is to be recreated, in a preliminary stage, intermediate data series, as
specified in `m.data_series`, are loaded in levels using `load_data_levels`. See
`?load_data_levels` for more details.

Then, the series in levels are transformed as specified in `m.data_transforms`. See
`?transform_data` for more details.

The resulting DataFrame is saved to disk as `data_<yymmdd>.csv` and returned to the caller.  
"""
function load_data(m::AbstractModel; try_disk::Bool = true, verbose::Symbol=:low)
    recreate_data = false

    # Check if already downloaded
    if try_disk && has_saved_data(m)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            print("Reading dataset from disk...")
        end
        df = read_data(m)
        if isvalid_data(m, df)
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                println("dataset from disk valid")
            end
        else
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                println("dataset from disk not valid")
            end
            recreate_data = true
        end
    else
        recreate_data = true
    end

    # Download routines
    if recreate_data
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Creating dataset...")
        end
        df = load_data_levels(m; verbose=verbose)
        df = transform_data(m, df; verbose=verbose)

        # Ensure that only appropriate rows make it into the returned DataFrame.
        start_date = date_presample_start(m)
        end_date   = date_zlb_end(m)
        df = df[start_date .<= df[:, :date] .<= end_date, :]

        save_data(m, df)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("dataset creation successful")
        end
    end

    return df
end

"""
```
load_data_levels(m::AbstractModel; verbose::Symbol=:low)
```

Load data in levels by appealing to the data sources specified for the model. Data from FRED
is loaded first, by default; then, merge other custom data sources.

Check on disk in `inpath(m, "data")` datasets, of the correct vintage, corresponding to the
ones in `keys(m.data_series)`. Load the appropriate data series (specified in
`m.data_series[source]`) for each data source. 
    
To accomodate growth rates and other similar transformations, more rows of data may be
downloaded than otherwise specified by the date model settings. (By the end of the process,
these rows will have been dropped.)

Data from FRED (i.e. the `:fred` data source) are treated separately. These are downloaded
using `load_fred_data`. See `?load_fred_data` for more details.

Data from non-FRED data sources are read from disk, verified, and merged.
"""
function load_data_levels(m::AbstractModel; verbose::Symbol=:low)
    # Start two quarters further back than `start_date` as we need these additional
    # quarters to compute differences.
    start_date = date_presample_start(m) - Dates.Month(6)
    end_date = date_zlb_end(m)

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

        # Skip FRED sources, which have already been handled
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
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                println("Reading data from $file...")
            end

            # Read in dataset and check that the file contains data for the proper dates
            addl_data = readtable(file)

            # Convert dates from strings to quarter-end dates for date arithmetic
            format_dates!(:date, addl_data)

            # Warn on sources with incomplete data; missing data will be replaced with NaN
            # during merge.
            if !in(lastdayofquarter(start_date), addl_data[:date]) ||
               !in(lastdayofquarter(end_date), addl_data[:date])
                warn("$file does not contain the entire date range specified; NaNs used.")
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

Save `df` to disk as CSV. File is located in `inpath(m, "data")`.
"""
function save_data(m::AbstractModel, df::DataFrame)
    vint = data_vintage(m)
    filename = inpath(m, "data", "data_$vint.csv")
    writetable(filename, df)
end

"""
```
has_saved_data(m::AbstractModel)
```

Determine if there is a saved dataset on disk for the required vintage.
"""
function has_saved_data(m::AbstractModel)
    vint = data_vintage(m)
    filename = inpath(m, "data", "data_$vint.csv")
    isfile(filename)
end

"""
```
read_data(m::AbstractModel)
```

Read CSV from disk as DataFrame. File is located in `inpath(m, "data")`.
"""
function read_data(m::AbstractModel)
    vint     = data_vintage(m)
    filename = inpath(m, "data", "data_$vint.csv")
    df       = readtable(filename)

    # Convert date column from string to Date
    df[:date] = map(Date, df[:date])

    return df
end

"""
```
isvalid_data(m::AbstractModel, df::DataFrame)
```

Return if dataset is valid for this model, ensuring that all observables are contained and
that all quarters between the beginning of the presample and the end of the mainsample are
contained.
"""
function isvalid_data(m::AbstractModel, df::DataFrame)
    valid = true

    # Ensure that every series in m_series is present in df_series
    m_series = collect(keys(m.data_transforms))
    df_series = names(df)
    coldiff = setdiff(m_series, df_series)
    valid = valid && isempty(coldiff)
    if !isempty(coldiff)
        println("Columns of 'df' do not match expected.")
        println(coldiff)
    end

    # Ensure the dates between date_presample_start and date_zlb_end are contained.
    actual_dates = df[:date]
    expected_dates = get_quarter_ends(date_presample_start(m), date_zlb_end(m))
    datesdiff = setdiff(expected_dates, actual_dates)
    valid = valid && isempty(datesdiff)
    if !isempty(datesdiff)
        println("Dates of 'df' do not match expected.")
        println(datesdiff)
    end

    return valid
end

"""
```
df_to_matrix(m::AbstractModel, df::DataFrame)
```

Return `df`, converted to matrix of floats, and discard date column. Also ensure data are
sorted by date and that rows outside of sample are discarded. The output of this function is
suitable for direct use in `estimate`, `posterior`, etc.
"""
function df_to_matrix(m::AbstractModel, df::DataFrame)
    # Sort rows by date and discard rows outside of sample
    df1 = sort(df; cols=[:date])
    t_start = find(df1[:date] .== date_presample_start(m))[1]
    t_end = find(df1[:date] .== date_zlb_end(m))[1]

    # Discard columns not used.
    #TODO sort columns as well according to observables
    datecol = df1.colindex[:date]
    colinds = setdiff(1:size(df1,2), datecol)

    df1 = df1[t_start:t_end, colinds]

    return convert(Matrix{Float64}, df1)
end
