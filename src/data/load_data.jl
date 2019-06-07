"""
```
load_data(m::AbstractModel; try_disk::Bool = true, verbose::Symbol = :low)
```

Create a DataFrame with all data series for this model, fully transformed.

First, check the disk to see if a valid dataset is already stored in `inpath(m, \"data\")`. A
dataset is valid if every series in `m.observable_mappings` is present and the entire sample is
contained (from `date_presample_start` to `date_mainsample_end`. If no valid dataset is
already stored, the dataset will be recreated. This check can be eliminated by passing
`try_disk=false`.

If the dataset is to be recreated, in a preliminary stage,
intermediate data series as specified in `m.observable_mappings` are
loaded in levels using `load_data_levels`. See `?load_data_levels` for
more details.

Then, the series in levels are transformed as specified in `m.observable_mappings`. See
`?transform_data` for more details.

If `m.testing` is false, then the resulting DataFrame is saved to disk as `data_<yymmdd>.csv`.
The data are then returned to the caller.
"""
function load_data(m::AbstractModel; cond_type::Symbol = :none, try_disk::Bool = true, verbose::Symbol=:low)
    recreate_data = false

    # Check if already downloaded
    if try_disk && has_saved_data(m; cond_type=cond_type)
        filename = get_data_filename(m, cond_type)
        print(verbose, :low, "Reading dataset $(filename) from disk...")
        df = read_data(m; cond_type = cond_type)
        if isvalid_data(m, df; cond_type = cond_type)
            println(verbose, :low, "dataset from disk valid")
        else
            println(verbose, :low, "dataset from disk not valid")
            recreate_data = true
        end
    else
        recreate_data = true
    end

    # Download routines
    if recreate_data
        println(verbose, :low, "Creating dataset...")

        levels = load_data_levels(m; verbose=verbose)
        if cond_type in [:semi, :full]
            cond_levels = load_cond_data_levels(m; verbose=verbose)
            levels, cond_levels = reconcile_column_names(levels, cond_levels)
            levels = vcat(levels, cond_levels)
        end
        df = transform_data(m, levels; cond_type=cond_type, verbose=verbose)

        # Ensure that only appropriate rows make it into the returned DataFrame.
        start_date = date_presample_start(m)
        end_date   = if cond_type in [:semi, :full]
            date_conditional_end(m)
        else
            date_mainsample_end(m)
        end
        df = df[start_date .<= df[:date] .<= end_date, :]

        if !m.testing
            save_data(m, df; cond_type=cond_type)
        end
        println(verbose, :low, "dataset creation successful")

        missing_cond_vars!(m, df; cond_type = cond_type)

        # check that dataset is valid
        isvalid_data(m, df)
    end

    return df
end

"""
```
load_data_levels(m::AbstractModel; verbose::Symbol=:low)
```

Load data in levels by appealing to the data sources specified for the model. Data from FRED
is loaded first, by default; then, merge other custom data sources.

Check on disk in `inpath(m, \"data\")` datasets, of the correct
vintage, corresponding to the ones required by the entries in
`m.observable_mappings`. Load the appropriate data series (specified
in `m.observable_mappings[key].input_series`) for each data source.

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
    end_date = date_mainsample_end(m)

    # Parse m.observable_mappings for data series
    data_series = parse_data_series(m)

    # Load FRED data
    df = load_fred_data(m; start_date=firstdayofquarter(start_date), end_date=end_date, verbose=verbose)


    # Set ois series to load
    if n_anticipated_shocks(m) > 0
        data_series[:OIS] = [Symbol("ant$i") for i in 1:n_anticipated_shocks(m)]
    end

    # For each additional source, search for the file with the proper name. Open
    # it, read it in, and merge it with fred_series

    vint = data_vintage(m)

    for source in keys(data_series)

        # Check that this source is actually used
        mnemonics = data_series[source]
        if isempty(mnemonics)
            @warn "No series were specified from $(string(source))"
            continue
        end

        # Skip FRED sources, which have already been handled
        # Conditional data are handled in `load_cond_data_levels`
        if source in [:FRED, :conditional]
            continue
        end

        # Read and merge data from this source
        file = inpath(m, "raw", "$(lowercase(string(source)))_$vint.csv")

        if isfile(file)
            println(verbose, :low, "Reading data from $file...")

            # Read in dataset and check that the file contains data for the proper dates
            addl_data = CSV.read(file, copycols=true)

            # Convert dates from strings to quarter-end dates for date arithmetic
            format_dates!(:date, addl_data)

            # Warn on sources with incomplete data; missing data will be replaced with missing
            # during merge.
            if !in(lastdayofquarter(start_date), addl_data[:date]) ||
                !in(lastdayofquarter(end_date), addl_data[:date])

                @warn "$file does not contain the entire date range specified; missings used."
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
            # If series not found, use all missings
            addl_data = DataFrame(fill(missing, (size(df,1), length(mnemonics))))
            names!(addl_data, mnemonics)
            df = hcat(df, addl_data)
            @warn "$file was not found; missings used"
        end
    end

    sort!(df, :date)

    # print population level data to a file
    if !m.testing
        filename = inpath(m, "raw", "population_data_levels_$vint.csv")
        mnemonic = parse_population_mnemonic(m)[1]
        if !isnull(mnemonic)
            CSV.write(filename, df[[:date, get(mnemonic)]])
        end
    end

    return df
end

"""
```
load_cond_data_levels(m::AbstractModel; verbose::Symbol=:low)
```

Check on disk in `inpath(m, \"cond\")` for a conditional dataset (in levels) of the correct
vintage and load it.

The following series are also loaded from `inpath(m, \"raw\")` and either
appended or merged into the conditional data:

- The last period of (unconditional) data in levels
  (`data_levels_<yymmdd>.csv`), used to calculate growth rates
- The first period of forecasted population
  (`population_forecast_<yymmdd>.csv`), used for per-capita calculations
"""
function load_cond_data_levels(m::AbstractModel; verbose::Symbol=:low)

    # Prepare file name
    cond_vint = cond_vintage(m)
    cond_idno = lpad(string(cond_id(m)), 2, string(0)) # print as 2 digits
    file = inpath(m, "cond", "cond_cdid=" * cond_idno * "_cdvt=" * cond_vint * ".csv")

    if isfile(file)
        println(verbose, :low, "Reading conditional data from $file...")

        # Read data
        cond_df = CSV.read(file, copycols=true)
        format_dates!(:date, cond_df)

        date_cond_end = cond_df[end, :date]

        # Use population forecast as population data
        population_forecast_file = inpath(m, "raw", "population_forecast_" * data_vintage(m) * ".csv")
        if isfile(population_forecast_file) && !isnull(get_setting(m, :population_mnemonic))
            pop_forecast = CSV.read(population_forecast_file, copycols=true)

            population_mnemonic = get(parse_population_mnemonic(m)[1])
            rename!(pop_forecast, :POPULATION =>  population_mnemonic)
            #DSGE.na2nan!(pop_forecast)
            DSGE.format_dates!(:date, pop_forecast)

            cond_df = join(cond_df, pop_forecast, on=:date, kind=:left)

            # Turn NAs into NaNs
            #na2nan!(cond_df)
            sort!(cond_df, :date)

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
    filename = get_data_filename(m, cond_type)
    CSV.write(filename, df)
end

"""
```
has_saved_data(m::AbstractModel; cond_type::Symbol = :none)
```

Determine if there is a saved dataset on disk for the required vintage and
conditional type.
"""
function has_saved_data(m::AbstractModel; cond_type::Symbol = :none)
    filename = get_data_filename(m, cond_type)
    return isfile(filename)
end

"""
```
read_data(m::AbstractModel; cond_type::Symbol = :none)
```

Read CSV from disk as DataFrame. File is located in `inpath(m, \"data\")`.
"""
function read_data(m::AbstractModel; cond_type::Symbol = :none)
    filename = get_data_filename(m, cond_type)
    df       = CSV.read(filename, copycols=true)

    # Convert date column from string to Date
    df[:date] = map(Date, df[:date])

    missing_cond_vars!(m, df; cond_type = cond_type)

    return df
end

"""
```
isvalid_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
```

Return if dataset is valid for this model, ensuring that all observables are contained and
that all quarters between the beginning of the presample and the end of the mainsample are
contained. Also checks to make sure that expected interest rate data is available if `n_anticipated_shocks(m) > 0`.
"""
function isvalid_data(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
    valid = true

    # Ensure that every series in m_series is present in df_series
    m_series = collect(keys(m.observable_mappings))
    df_series = names(df)
    coldiff = setdiff(m_series, df_series)
    valid = valid && isempty(coldiff)
    if !isempty(coldiff)
        println("Columns of 'df' do not match expected.")
        println(coldiff)
    end

    # Ensure the dates between date_presample_start and date_mainsample_end are contained.
    actual_dates = df[:date]

    start_date = date_presample_start(m)
    end_date   = if cond_type in [:semi, :full]
        date_conditional_end(m)
    else
        date_mainsample_end(m)
    end
    expected_dates = get_quarter_ends(start_date, end_date)
    datesdiff = setdiff(expected_dates, actual_dates)

    valid = valid && isempty(datesdiff)
    if !isempty(datesdiff)
        println("Dates of 'df' do not match expected.")
        println(datesdiff)
    end

    # Ensure that no series is all missing
    for col in setdiff(names(df), [:date])
        if all(ismissing.(df[col]))
            @warn "df[$col] is all missing."
        end
    end

    return valid
end

"""
```
df_to_matrix(m, df; cond_type = :none, in_sample = true)
```

Return `df`, converted to matrix of floats, and discard date column. Also ensure
that rows are sorted by date and columns by `m.observables`, with the option to specify
whether or not the out of sample rows are discarded. The output of this
function is suitable for direct use in `estimate`, `posterior`, etc.

### Keyword Arguments:
- `include_presample::Bool`: indicates whether or not there are presample periods.
- `in_sample::Bool`: indicates whether or not to discard rows that are out of sample. Set this flag to false in
the case that you are calling filter_shocks! in the scenarios codebase.
"""
function df_to_matrix(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none, include_presample::Bool = true, in_sample::Bool = true)
    # Sort rows by date
    df1 = sort(df, :date)

    if in_sample
        start_date = if include_presample
            date_presample_start(m)
        else
            date_mainsample_start(m)
        end
        end_date = if cond_type in [:semi, :full]
            date_conditional_end(m)
        else
            date_mainsample_end(m)
        end
        df1 = df1[start_date .<= df[:date] .<= end_date, :]
    end

    # Discard columns not used
    cols = collect(keys(m.observables))
    sort!(cols, by = x -> m.observables[x])
    df1 = df1[cols]

    return permutedims(Float64.(collect(Missings.replace(convert(Matrix{Union{Missing, Float64}}, df1), NaN))))
end

"""
```
data_to_df(m, data, start_date)
```

Create a `DataFrame` out of the matrix `data`, including a `:date` column
beginning in `start_date`.  Variable names and indices are obtained from
`m.observables`.
"""
function data_to_df(m::AbstractModel, data::Matrix{T}, start_date::Date) where T<:AbstractFloat
    # Check number of rows = number of observables
    nobs = n_observables(m)
    @assert size(data, 1) == nobs "Number of rows of data matrix ($(size(data, 1))) must equal number of observables ($nobs)"

    # Initialize DataFrame and add dates
    nperiods = size(data, 2)
    end_date = iterate_quarters(start_date, nperiods - 1)
    dates = quarter_range(start_date, end_date)
    df = DataFrame(date = dates)

    # Add observables
    for var in keys(m.observables)
        ind = m.observables[var]
        df[var] = vec(data[ind, :])
    end

    return df
end

"""
```
parse_data_series(m::AbstractModel)
```

Parse `m.observable_mappings` for the data sources and mnemonics to read in.

Returns a `Dict{Symbol, Vector{Symbol}}` mapping sources => mnemonics found in that data file.
"""
function parse_data_series(m::AbstractModel)

    data_series = Dict{Symbol, Vector{Symbol}}()

    # Parse vector of observable mappings into data_series dictionary
    for obs in values(m.observable_mappings)
        for series in obs.input_series
            mnemonic, source = map(Symbol, split(string(series), DSGE_DATASERIES_DELIM))

            if !in(source, keys(data_series))
                data_series[source] = Vector{Symbol}()
            end

            if !in(mnemonic, data_series[source])
                push!(data_series[source], mnemonic)
            end
        end
    end
    data_series
end

"""
```
read_population_data(m; verbose = :low)

read_population_data(filename; verbose = :low)
```

Read in population data stored in levels, either from
`inpath(m, \"raw\", \"population_data_levels_[vint].csv\"`) or `filename`.
"""
function read_population_data(m::AbstractModel; verbose::Symbol = :low)
    vint = data_vintage(m)
    filename = inpath(m, "raw", "population_data_levels_" * vint * ".csv")
    read_population_data(filename; verbose = verbose)
end

function read_population_data(filename::String; verbose::Symbol = :low)
    println(verbose, :low, "Reading population data from $filename...")

    df = CSV.read(filename, copycols=true)
    DSGE.format_dates!(:date, df)
    sort!(df, :date)

    return df
end

"""
```
read_population_forecast(m; verbose = :low)

read_population_forecast(filename, population_mnemonic, last_recorded_date; verbose = :low)
```

Read in population forecast in levels, either from
`inpath(m, \"raw\", \"population_forecast_[vint].csv\")` or `filename`.
If that file does not exist, return an empty `DataFrame`.

"""
function read_population_forecast(m::AbstractModel; verbose::Symbol = :low)
    population_forecast_file = inpath(m, "raw", "population_forecast_" * data_vintage(m) * ".csv")
    population_mnemonic = parse_population_mnemonic(m)[1]

    if isnull(population_mnemonic)
        error("No population mnemonic provided")
    else
        read_population_forecast(population_forecast_file, get(population_mnemonic); verbose = verbose)
    end
end

function read_population_forecast(filename::String, population_mnemonic::Symbol;
                                  verbose::Symbol = :low)
    if isfile(filename)
        println(verbose, :low, "Loading population forecast from $filename...")

        df = CSV.read(filename, copycols=true)
        rename!(df, :POPULATION => population_mnemonic)
        DSGE.format_dates!(:date, df)
        sort!(df, :date)

        return df[[:date, population_mnemonic]]
    else
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @warn "No population forecast found"
        end
        return DataFrame()
    end
end
