"""
```
load_data(m::AbstractDSGEModel; try_disk::Bool = true, verbose::Symbol = :low,
          check_empty_columns::Bool = true, summary_statistics::Symbol = :low)
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

The keyword `check_empty_columns` throws an error whenever a column is completely empty
in the loaded data set if it is set to true.

The keyword `summary_statistics` prints out a variety of summary statistics
on the loaded data. When set to :low, we print only the number of
missing/NaNs for each data series. When set to :high, we also print
means, standard deviations,
"""
function load_data(m::AbstractDSGEModel; cond_type::Symbol = :none, try_disk::Bool = true,
                   verbose::Symbol=:low, check_empty_columns::Bool = true,
                   summary_statistics::Symbol = :low, add_vals = (false, Date(2020,12,31)))
    recreate_data = false

    # Check if already downloaded
    if try_disk && has_saved_data(m; cond_type=cond_type)
        filename = get_data_filename(m, cond_type)
        print(verbose, :low, "Reading dataset $(filename) from disk...")
        df = read_data(m; cond_type = cond_type, check_empty_columns = check_empty_columns)
        if isvalid_data(m, df; cond_type = cond_type, check_empty_columns = check_empty_columns)
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

        levels = load_data_levels(m; verbose=verbose, add_vals = add_vals)
        if cond_type in [:semi, :full]
            cond_levels = load_cond_data_levels(m; verbose=verbose)
            levels, cond_levels = reconcile_column_names(levels, cond_levels)
            levels = vcat(levels, cond_levels)
        end
        df = transform_data(m, levels; cond_type=cond_type, verbose=verbose)

        if :obs_nominalrate1 in cond_semi_names(m) || :obs_nominalrate1 in cond_full_names(m)
            ois_data = CSV.read(inpath(m, "raw", "ois_$(data_vintage(m)).csv"), DataFrame, copycols = true)
            dates = DSGE.get_quarter_ends(iterate_quarters(date_mainsample_end(m), 1), date_conditional_end(m))
            n_cond = length(dates)

            ois_data_want = ois_data[date_mainsample_end(m) .< ois_data[!, :date] .<= date_conditional_end(m), [Symbol("ant$i") for i in 1:n_mon_anticipated_shocks(m)]]

            df[date_mainsample_end(m) .< df[!, :date] .<= date_conditional_end(m), [Symbol("obs_nominalrate$i") for i in 1:n_mon_anticipated_shocks(m)]] .= Matrix{Float64}(ois_data_want)
        end

        # Ensure that only appropriate rows make it into the returned DataFrame.
        start_date = date_presample_start(m)
        end_date   = if cond_type in [:semi, :full]
            date_conditional_end(m)
        else
            date_mainsample_end(m)
        end
        df = df[start_date .<= df[!,:date] .<= end_date, :]

        missing_cond_vars!(m, df; cond_type = cond_type,
                           check_empty_columns = check_empty_columns)

        # check that dataset is valid
        isvalid_data(m, df; check_empty_columns = check_empty_columns)

        if !m.testing
            save_data(m, df; cond_type=cond_type)
        end
        println(verbose, :low, "dataset creation successful")

        # print summary statistics
        if summary_statistics == :low || summary_statistics == :high
            str_nondate_names = [string(name) for name in propertynames(df[:,2:end])]
            freq_nan_empty = zeros(size(df,2) - 1)
            for (colnum, name) in enumerate(propertynames(df[:,2:end]))
                is_missing_in_col = ismissing.(df[!,name])
                is_nan_in_col = isnan.(df[!,name][.!is_missing_in_col])
                n_miss = count(is_missing_in_col)
                n_nan  = count(is_nan_in_col)
                freq_nan_empty[colnum] = (n_nan + n_miss) / length(is_missing_in_col)
                println("$(name), Frequency of missing/NaNs: $(freq_nan_empty[colnum])")
                if summary_statistics == :high
                    colmean = mean(df[!,name][.!is_missing_in_col][.!is_nan_in_col])
                    colstd  = std(df[!,name][.!is_missing_in_col][.!is_nan_in_col])
                    println("$(name), Column Mean: $(colmean)")
                    println("$(name), Standard deviation: $(colstd)")
                end
            end
        end
    end

    return df
end

"""
```
load_data_levels(m::AbstractDSGEModel; verbose::Symbol=:low)
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
function load_data_levels(m::AbstractDSGEModel; verbose::Symbol=:low,
                          add_vals = (false,Date(2020,12,31)))
    # Start two quarters further back than `start_date` as we need these additional
    # quarters to compute differences.
    start_date = date_presample_start(m) - Dates.Month(6)
    end_date = date_mainsample_end(m)

    # Parse m.observable_mappings for data series
    data_series = parse_data_series(m)

    # Load FRED data
    df = load_fred_data(m; start_date=firstdayofquarter(start_date),
                        end_date=end_date, verbose=verbose)

    if add_vals[1]
        temp2 = CSV.read(get_setting(m, :cond_filename), DataFrame)
        for k in names(temp2)[2:end]
            if k in names(df)
                df[df[:date] .== add_vals[2],k] = temp2[1,k]
            end
        end
    end

    # Set ois series to load
    if n_mon_anticipated_shocks(m) > 0
        if get_setting(m, :rate_expectations_source) == :ois
            data_series[:OIS] = [Symbol("ant$i") for i in 1:n_mon_anticipated_shocks(m)]
        elseif get_setting(m, :rate_expectations_source) == :bluechip
            data_series[:bluechip] = [Symbol("ant$i") for i in 1:n_mon_anticipated_shocks(m)]
        end
    end
    try
        if n_z_anticipated_shocks(m) > 0
            data_series[:Z] = [Symbol("z$i") for i in 1:n_z_anticipated_shocks(m)]
        end
    catch err
        if isa(err, KeyError)
            nothing
        else
            rethrow(err)
        end
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
            addl_data = CSV.read(file, DataFrame, copycols = true)

            # Convert dates from strings to quarter-end dates for date arithmetic
            format_dates!(:date, addl_data)

            # Warn on sources with incomplete data; missing data will be replaced with missing
            # during merge.
            if !in(lastdayofquarter(start_date), addl_data[!,:date]) ||
                !in(lastdayofquarter(end_date), addl_data[!,:date])

                @warn "$file does not contain the entire date range specified; missings used."
            end

            # Make sure each mnemonic that was specified is present
            for series in mnemonics
                if !in(series, propertynames(addl_data))
                    error("$(string(series)) is missing from $file.")
                end
            end

            # Extract just the columns and rows of the dataset we want, and merge them with
            # data
            cols = [:date; mnemonics]
            rows = start_date .<= addl_data[!,:date] .<= end_date

            addl_data = addl_data[rows, cols]
            if isdefined(DataFrames, :outerjoin)
                df = outerjoin(df, addl_data, on=:date)
            else
                df = join(df, addl_data, on=:date, kind = :outer)
            end
        else
            # If series not found, use all missings
            addl_data = DataFrame(fill(missing, (size(df,1), length(mnemonics))))
            if isdefined(DataFrames, :rename!)
                rename!(addl_data, mnemonics)
            else
                names!(addl_data, mnemonics)
            end
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
            CSV.write(filename, df[!,[:date, get(mnemonic)]])
        end
    end
    return df
end

"""
```
load_cond_data_levels(m::AbstractDSGEModel; verbose::Symbol=:low)
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
function load_cond_data_levels(m::AbstractDSGEModel; verbose::Symbol=:low)

    # Prepare file name
    cond_vint = cond_vintage(m)
    cond_idno = lpad(string(cond_id(m)), 2, string(0)) # print as 2 digits
    file = inpath(m, "cond", "cond_cdid=" * cond_idno * "_cdvt=" * cond_vint * ".csv")

    if isfile(file)
        println(verbose, :low, "Reading conditional data from $file...")

        # Read data
        cond_df = CSV.read(file, DataFrame, copycols = true)
        format_dates!(:date, cond_df)

        date_cond_end = cond_df[end, :date]

        # Use population forecast as population data
        population_forecast_file = inpath(m, "raw", "population_forecast_" * data_vintage(m) * ".csv")
        if isfile(population_forecast_file) && !isnull(get_setting(m, :population_mnemonic))
            pop_forecast = CSV.read(population_forecast_file, DataFrame, copycols = true)

            population_mnemonic = get(parse_population_mnemonic(m)[1])
            rename!(pop_forecast, :POPULATION =>  population_mnemonic)
            # na2nan!(pop_forecast) # Removed b/c DataFrames uses missing instead of NA, and missings are already handled
            format_dates!(:date, pop_forecast)

            cond_df = if isdefined(DataFrames, :leftjoin) # left joins using `join` is deprecated in DataFrames v0.21 (and higher)
                leftjoin(cond_df, pop_forecast, on = :date)
            else
                join(cond_df, pop_forecast, on = :date, kind = :left)
            end

            # Turn NAs into NaNs
            # na2nan!(cond_df) # Removed b/c DataFrames uses missing instead of NA, and missings are already handled

            # Make sure the data is ordered by the date
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
save_data(m::AbstractDSGEModel, df::DataFrame; cond_type::Symbol = :none)
```

Save `df` to disk as CSV. File is located in `inpath(m, \"data\")`.
"""
function save_data(m::AbstractDSGEModel, df::DataFrame; cond_type::Symbol = :none)
    filename = get_data_filename(m, cond_type)
    CSV.write(filename, df)
end

"""
```
has_saved_data(m::AbstractDSGEModel; cond_type::Symbol = :none)
```

Determine if there is a saved dataset on disk for the required vintage and
conditional type.
"""
function has_saved_data(m::AbstractDSGEModel; cond_type::Symbol = :none)
    filename = get_data_filename(m, cond_type)
    return isfile(filename)
end

"""
```
read_data(m::AbstractDSGEModel; cond_type::Symbol = :none)
```

Read CSV from disk as DataFrame. File is located in `inpath(m, \"data\")`.
"""
function read_data(m::AbstractDSGEModel; cond_type::Symbol = :none, check_empty_columns::Bool = true)
    filename = get_data_filename(m, cond_type)
    df       = CSV.read(filename, DataFrame, copycols = true)

    # Convert date column from string to Date
    df[!,:date] = map(Date, df[!,:date])

    missing_cond_vars!(m, df; cond_type = cond_type, check_empty_columns = check_empty_columns)

    return df
end

"""
```
isvalid_data(m::AbstractDSGEModel, df::DataFrame; cond_type::Symbol = :none,
    check_empty_columns::Bool = true)
```

Return if dataset is valid for this model, ensuring that all observables are contained and
that all quarters between the beginning of the presample and the end of the mainsample are
contained. Also checks to make sure that expected interest rate data is available if `n_mon_anticipated_shocks(m) > 0`.
"""
function isvalid_data(m::AbstractDSGEModel, df::DataFrame; cond_type::Symbol = :none,
                      check_empty_columns::Bool = true)
    valid = true

    # Ensure that every series in m_series is present in df_series
    m_series = collect(keys(m.observable_mappings))
    df_series = propertynames(df)
    coldiff = setdiff(m_series, df_series)
    valid = valid && isempty(coldiff)
    if !isempty(coldiff)
        println("Columns of 'df' do not match expected.")
        println(coldiff)
    end

    # Ensure the dates between date_presample_start and date_mainsample_end are contained.
    actual_dates = df[!,:date]

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

    # Ensure that no series is all missing or NaNs
    if check_empty_columns
        empty_cols = Vector{String}(undef,0)
        for name in setdiff(propertynames(df), [:date])
            is_missing_in_col = ismissing.(df[!, name])
            is_nan_in_col = isnan.(df[!, name][.!is_missing_in_col])
            if sum(vcat(is_missing_in_col, is_nan_in_col)) == length(df[!, name])
                push!(empty_cols, string(name) * ", ")
            end
        end
        as_str = join(empty_cols)
        if !isempty(empty_cols)
            error("Column(s) $(as_str[1:end-2]) have only NaNs and/or missings.")
        end
    else
        for col in setdiff(propertynames(df), [:date])
            if all(ismissing.(df[!,col]))
                @warn "df[$col] is all missing."
            else
                not_missing_inds = .!(ismissing.(df[!, col]))
                if all(isnan.(df[not_missing_inds, col]))
                    @warn "df[$col] is entirely NaNs and/or missings."
                end
            end
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
function df_to_matrix(m::Union{AbstractDSGEModel,AbstractVARModel}, df::DataFrame; cond_type::Symbol = :none, include_presample::Bool = true, in_sample::Bool = true)
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
        df1 = df1[start_date .<= df[!,:date] .<= end_date, :]
    end

    # Discard columns not used
    cols = collect(keys(get_observables(m)))
    sort!(cols, by = x -> get_observables(m)[x])
    df1 = df1[!,cols]

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
function data_to_df(m::AbstractDSGEModel, data::Matrix{T}, start_date::Date) where T<:AbstractFloat
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
        df[!, var] = vec(data[ind, :])
    end

    return df
end

"""
```
parse_data_series(m::AbstractDSGEModel)
```

Parse `m.observable_mappings` for the data sources and mnemonics to read in.

Returns a `Dict{Symbol, Vector{Symbol}}` mapping sources => mnemonics found in that data file.
"""
function parse_data_series(m::AbstractDSGEModel)

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
function read_population_data(m::AbstractDSGEModel; verbose::Symbol = :low)
    vint = data_vintage(m)
    filename = inpath(m, "raw", "population_data_levels_" * vint * ".csv")
    read_population_data(filename; verbose = verbose)
end

function read_population_data(filename::String; verbose::Symbol = :low)
    println(verbose, :low, "Reading population data from $filename...")

    df = CSV.read(filename, DataFrame, copycols = true)
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
function read_population_forecast(m::AbstractDSGEModel; verbose::Symbol = :low)
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

        df = CSV.read(filename, DataFrame, copycols = true)
        rename!(df, :POPULATION => population_mnemonic)
        DSGE.format_dates!(:date, df)
        sort!(df, :date)

        return df[!,[:date, population_mnemonic]]
    else
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @warn "No population forecast found"
        end
        return DataFrame()
    end
end

# Takes  individual forecast
function construct_bluechip_data(m::AbstractModel, raw_forecasts::String, num_anticipated_shocks::Int, save_bluechip_rate_expectations_here::String)
    df = DataFrame(date = DSGE.quarter_range(date_zlb_start(m), date_zlb_end(m)))
    for h = 1:num_anticipated_shocks
        df[!, Symbol("ant$h")] = NaN
    end
    for date in DSGE.quarter_range(date_zlb_start(m), date_zlb_end(m))
        y, q = Dates.year(date), Dates.quarterofyear(date)
        bluechip = load_reference_forecast(m, y, q, raw_forecasts, num_anticipated_shocks)
        df_date_ind = findfirst(df[!, :date], date)
        for h=1:num_anticipated_shocks
            df[df_date_ind, Symbol("ant$h")] =
                # We use quarterly interest rates but bluechip forecasts annualized
                bluechip[h]/4
        end
    end
    CSV.write(save_bluechip_rate_expectations_here, df)
    print("Saved $(save_bluechip_rate_expectations_here)")
    return df
end

function load_reference_forecast(m::AbstractModel, year::Int, quarter::Int, file::String, num_anticipated_shocks::Int)
    # Dates
    start_date = DSGE.quartertodate(string(year)*"q"*string(quarter))
    end_date   = DSGE.iterate_quarters(start_date, num_anticipated_shocks)
    dates      = DSGE.quarter_range(start_date, end_date)

    df = DataFrame(date = dates)

    series = CSV.read(file, DataFrame)

    datestr = reference_forecast_vintage(year, quarter, :bluechip)
    date_index = something(findfirst(isequal(datestr), series[:date]), 0)

    horizon = num_anticipated_shocks
    # Assign t-quarters-ahead forecast of var to df[t, var]
    forecast = NaN*zeros(horizon)

    horizon = minimum([horizon, subtract_quarters(get_setting(m, :date_zlb_end), DSGE.quartertodate(string(year)*"q"*string(quarter))) + 1])
    for t = 1:horizon
        forecast[t] = series[date_index, Symbol("bluechip_nominalrate_", t)]
    end

    # Assign forecast of var to column of df
    #df[:obs_nominalrate] = forecast

    return forecast
end


function reference_forecast_vintage(year::Int, quarter::Int, reference_forecast::Symbol)
    # We compare to Jaunary, April, July, and October Blue Chip forecasts (bluechip_forecast_month = 1 in Realtime code
        bluechip_forecast_month = 1
        release_quarter = quarter % 4 + 1
        release_year = if release_quarter == 1
            year + 1
        else
            year
        end
        release_month = 3*(release_quarter - 1) + bluechip_forecast_month
        release_day   = 10

    return Date(release_year, release_month, release_day)
end
