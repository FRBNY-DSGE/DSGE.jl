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

If `m.testing` is false, then the resulting DataFrame is saved to disk as `data_<yymmdd>.csv`.
The data are then returned to the caller.
"""
function load_data(m::AbstractModel; cond_type::Symbol = :none, try_disk::Bool = true, verbose::Symbol=:low)
    recreate_data = false

    # Check if already downloaded
    if try_disk && has_saved_data(m; cond_type=cond_type)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            print("Reading dataset from disk...")
        end
        df = read_data(m; cond_type = cond_type)
        if isvalid_data(m, df; cond_type = cond_type)
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

        levels = load_data_levels(m; verbose=verbose)
        if cond_type in [:semi, :full]
            cond_levels = load_cond_data_levels(m; verbose=verbose)
            levels = vcat(levels, cond_levels)
            na2nan!(levels)
        end
        df = transform_data(m, levels; cond_type=cond_type, verbose=verbose)

        # Ensure that only appropriate rows make it into the returned DataFrame.
        start_date = date_presample_start(m)
        end_date   = if cond_type in [:semi, :full]
            date_conditional_end(m)
        else
            date_zlb_end(m)
        end
        df = df[start_date .<= df[:, :date] .<= end_date, :]

        if !m.testing
            save_data(m, df; cond_type=cond_type)
        end
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

Check on disk in `inpath(m, \"data\")` datasets, of the correct vintage, corresponding to the
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

        # Check that this source is actually used
        mnemonics = m.data_series[source]
        if isempty(mnemonics)
            warn("No series were specified from $(string(source))")
            continue
        end

        # Skip FRED sources, which have already been handled
        # Conditional data are handled in `load_cond_data_levels`
        if source == :fred || source == :conditional
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
            warn("$file was not found; NaNs used.")
        end
    end
    
    # turn NAs into NaNs
    na2nan!(df)

    sort!(df, cols = :date)
end

"""
```
load_cond_data_levels(m::AbstractModel; verbose::Symbol=:low)
```

Load conditional data in levels.

Check on disk in `inpath(m, \"cond\")` for a conditional dataset of the correct
vintage. Load the appropriate data series, specified in
`m.data_series[:conditional]`.
    
The following series are also loaded from `inpath(m, \"data\")` and either
appended or merged into the conditional data:

- The last period of (unconditional) data in levels
  (`data_levels_<yymmdd>.csv`), used to calculate growth rates
- The first period of forecasted population
  (`population_forecast_<yymmdd>.csv`), used for per-capita calculations
"""
function load_cond_data_levels(m::AbstractModel; verbose::Symbol=:low)

    # Prepare file name
    cond_vint = get_setting(m, :cond_vintage)
    cond_id = get_setting(m, :cond_id)
    file = inpath(m, "cond", "cond_vint=$(cond_vint)_cdid=$(cond_id).csv")

    if isfile(file)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Reading conditional data from $file...")
        end

        # Read data
        cond_df = readtable(file)
        format_dates!(:date, cond_df)

        date_cond_end = cond_df[end, :date]

        # Make sure each mnemonic that was specified is present
        mnemonics = m.data_series[:conditional]
        for series in mnemonics
            if !in(series, names(cond_df))
                error("$(string(series)) is missing from $file.")
            end
        end

        # Use population forecast as population data
        population_forecast_file = inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
        if isfile(population_forecast_file)
            pop_forecast = readtable(population_forecast_file)

            population_mnemonic = get_setting(m, :population_mnemonic)
            rename!(pop_forecast, :POPULATION,  population_mnemonic)
            DSGE.na2nan!(pop_forecast)
            DSGE.format_dates!(:date, pop_forecast)

            cond_df = join(cond_df, pop_forecast, on=:date, kind=:left)

            # Turn NAs into NaNs
            na2nan!(cond_df)
            sort!(cond_df, cols = :date)

            return cond_df
        else
            error("Population forecast data in $population_forecast_file not found, but required to load conditional data")
        end
    else
        # If series not found, throw an error
        error("Conditional data in $file not found")
    end
end

"""
```
save_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
```

Save `df` to disk as CSV. File is located in `inpath(m, \"data\")`.
"""
function save_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
    vint = data_vintage(m)
    filestring = "data"
    if cond_type in [:semi, :full]
        filestring = filestring * "_cond=$cond_type"
    end
    filename = inpath(m, "data", "$(filestring)_$vint.csv")
    writetable(filename, df)
end

"""
```
has_saved_data(m::AbstractModel; cond_type::Symbol = :none)
```

Determine if there is a saved dataset on disk for the required vintage and
conditional type.
"""
function has_saved_data(m::AbstractModel; cond_type::Symbol = :none)
    vint = data_vintage(m)
    filestring = "data"
    if cond_type in [:semi, :full]
        filestring = filestring * "_cond=$cond_type"
    end
    filename = inpath(m, "data", "$(filestring)_$vint.csv")
    isfile(filename)
end

"""
```
read_data(m::AbstractModel; cond_type::Symbol = :none)
```

Read CSV from disk as DataFrame. File is located in `inpath(m, \"data\")`.
"""
function read_data(m::AbstractModel; cond_type::Symbol = :none)
    vint = data_vintage(m)
    filestring = "data"
    if cond_type in [:semi, :full]
        filestring = filestring * "_cond=$cond_type"
    end
    filename = inpath(m, "data", "$(filestring)_$vint.csv")
    df       = readtable(filename)

    # Convert date column from string to Date
    df[:date] = map(Date, df[:date])

    return df
end

"""
```
isvalid_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
```

Return if dataset is valid for this model, ensuring that all observables are contained and
that all quarters between the beginning of the presample and the end of the mainsample are
contained.
"""
function isvalid_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
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

    start_date = date_presample_start(m)
    end_date   = if cond_type in [:semi, :full]
        date_conditional_end(m)
    else
        date_zlb_end(m)
    end
    expected_dates = get_quarter_ends(start_date, end_date)
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
df_to_matrix(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
```

Return `df`, converted to matrix of floats, and discard date column. Also ensure data are
sorted by date and that rows outside of sample are discarded. The output of this function is
suitable for direct use in `estimate`, `posterior`, etc.
"""
function df_to_matrix(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
    # Sort rows by date and discard rows outside of sample
    df1 = sort(df; cols=[:date])

    start_date = date_presample_start(m)
    end_date   = if cond_type in [:semi, :full]
        date_conditional_end(m)
    else
        date_zlb_end(m)
    end
    df1 = df1[start_date .<= df1[:, :date] .<= end_date, :]

    # Discard columns not used.
    cols = collect(keys(m.observables))
    sort!(cols, by = x -> m.observables[x])
    df1 = df1[cols]

    return convert(Matrix{Float64}, df1)'
end
