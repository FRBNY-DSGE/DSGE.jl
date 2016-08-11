"""
```
load_cond_data(m::AbstractModel, cond_type::Symbol; verbose::Symbol = :low)
```

Create a DataFrame with conditional data for this model, fully transformed. The
`cond_type` can be one of:

- `:full`: uses forecasts of some economic observables (e.g. `:obs_gdp` and
  `:obs_corepce` in Model 990) and assumes that the averages of financial
  observables (`:obs_spread` and `:obs_nominalrate` in Model 990) observed in
  the days since the end of the last data period represent the observation for
  the next quarter. The observables used in conditional forecasts are specified
  in `cond_full_names(m)`.
- `:semi`: includes only financial data observed so far in the quarter
  (`:obs_spread` and `:obs_nominalrate` in Model 990). The observables used in
  conditional forecasts are specified in `cond_semi_names(m)`.
- `:none`: returns an empty DataFrame.

As in `load_data`, intermediate data series, as specified in
`m.data_series[:conditional]`, are loaded in levels using
`load_cond_data_levels`. See `?load_cond_data_levels` for more details.

Then, the series in levels are transformed as specified in `m.data_transforms`.
See `?transform_data` for more details.
"""
function load_cond_data(m::AbstractModel, cond_type::Symbol; verbose::Symbol = :low)
    
    if cond_type == :none
        return DataFrame()
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Creating conditional dataset...")
    end
    cond_levels, date_conditional_end = load_cond_data_levels(m; verbose=verbose)
    cond_df = transform_data(m, cond_levels; cond_type=cond_type, verbose=verbose)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("conditional dataset creation successful")
    end

    # Ensure that only appropriate rows make it into the returned DataFrame.
    start_date = date_forecast_start(m)
    end_date   = date_conditional_end

    cond_df = cond_df[start_date .<= cond_df[:, :date] .<= end_date, :]
    return cond_df
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
        
        # Need last period of (unconditional) data in levels to calculate growth rates
        vint = get_setting(m, :data_vintage)
        uncond_file = inpath(m, "data", "data_levels_$vint.csv")
        if !isfile(uncond_file)
            load_data_levels(m; verbose=verbose)
        end
        uncond_df = readtable(uncond_file)

        population_mnemonic = get_setting(m, :population_mnemonic)
        delete!(uncond_df, population_mnemonic)
        format_dates!(:date, uncond_df)

        uncond_lastperiod = uncond_df[uncond_df[:, :date] .== date_zlb_end(m), :]
        cond_df = vcat(uncond_lastperiod, cond_df)

        date_conditional_end = cond_df[end, :date]

        # Make sure each mnemonic that was specified is present
        mnemonics = m.data_series[:conditional]
        for series in mnemonics
            if !in(series, names(cond_df))
                error("$(string(series)) is missing from $file.")
            end
        end

        # Use population forecast as population data
        population_forecast_file = inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")

        pop_forecast = readtable(population_forecast_file)
        rename!(pop_forecast, :POPULATION,  population_mnemonic)
        DSGE.na2nan!(pop_forecast)
        DSGE.format_dates!(:date, pop_forecast)

        cond_df = join(cond_df, pop_forecast, on=:date, kind=:left)

        # turn NAs into NaNs
        na2nan!(cond_df)

        sort!(cond_df, cols = :date)

        return cond_df, date_conditional_end
    else
        # If series not found, throw an error
        error("Conditional data in $file not found")
    end
end