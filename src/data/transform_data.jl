"""
```
transform_data(m::AbstractModel, levels::DataFrame; cond_type::Symbol = :none,
    verbose::Symbol = :low)
```

Transform data loaded in levels and order columns appropriately for the DSGE
model. Returns DataFrame of transformed data.

The DataFrame `levels` is output from `load_data_levels`. The series in levels are
transformed as specified in `m.observable_mappings`.

- To prepare for per-capita transformations, population data are filtered using
  `hpfilter`. The series in `levels` to use as the population series is given by
  the `population_mnemonic` setting. If `use_population_forecast(m)`, a
  population forecast is appended to the recorded population levels before the
  filtering. Both filtered and unfiltered population levels and growth rates are
  added to the `levels` data frame.
- The transformations are applied for each series using the `levels` DataFrame
  as input.

Conditional data (identified by `cond_type in [:semi, :full]`) are handled
slightly differently: If `use_population_forecast(m)`, we drop the first period
of the population forecast because we treat the first forecast period
`date_forecast_start(m)` as if it were data. We also only apply transformations
for the observables given in `cond_full_names(m)` or `cond_semi_names(m)`.
"""
function transform_data(m::AbstractModel, levels::DataFrame; cond_type::Symbol = :none, verbose::Symbol = :low)

    n_obs, _ = size(levels)

    # Step 1: HP filter population forecasts, if they're being used
    population_mnemonic = parse_population_mnemonic(m)[1]
    if !isnull(population_mnemonic)
        population_forecast_levels = read_population_forecast(m; verbose = verbose)
        population_data, _ = transform_population_data(levels, population_forecast_levels,
                                                       get(population_mnemonic); verbose = verbose)

        levels = join(levels, population_data, on = :date, kind = :left)
        rename!(levels, [:filtered_population_recorded, :dlfiltered_population_recorded, :dlpopulation_recorded],
                [:filtered_population, :filtered_population_growth, :unfiltered_population_growth])
    end

    # Step 2: apply transformations to each series
    transformed = DataFrame()
    transformed[:date] = levels[:date]

    data_transforms = collect_data_transforms(m)

    for series in keys(data_transforms)
        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Transforming series $series...")
        end
        f = data_transforms[series]
        transformed[series] = f(levels)
    end

    sort!(transformed, cols = :date)

    return transformed
end

function collect_data_transforms(m; direction=:fwd)

    data_transforms = OrderedDict{Symbol,Function}()

    # Parse vector of observable mappings into data_transforms dictionary
    for obs in keys(m.observable_mappings)
        data_transforms[obs] = getfield(m.observable_mappings[obs], Symbol(string(direction) * "_transform"))
    end

    data_transforms
end

"""
```
transform_population_data(population_data, population_forecast,
    population_mnemonic; verbose = :low)
```

Load, HP-filter, and compute growth rates from population data in
levels. Optionally do the same for forecasts.

### Inputs

- `population_data`: pre-loaded DataFrame of historical population data
  containing the columns `:date` and `population_mnemonic`. Assumes this is
  sorted by date.
- `population_forecast`: pre-loaded `DataFrame` of population forecast
  containing the columns `:date` and `population_mnemonic`
- `population_mnemonic`: column name for population series in `population_data`
  and `population_forecast`

### Keyword Arguments

- `verbose`: one of `:none`, `:low`, or `:high`

### Output

A dictionary containing the following keys:

- `:filtered_population_recorded`: HP-filtered historical population series (levels)
- `:dlfiltered_population_recorded`: HP-filtered historical population series (growth rates)
- `:dlpopulation_recorded`: Non-filtered historical population series (growth rates)
- `:filtered_population_forecast`: HP-filtered population forecast series (levels)
- `:dlfiltered_population_forecast`: HP-filtered population forecast series (growth rates)
- `:dlpopulation_forecast`: Non-filtered historical population series (growth rates)

Note: the r\"*forecast\" fields will be empty if population_forecast_file is not provided.
"""
function transform_population_data(population_data::DataFrame, population_forecast::DataFrame,
                                   population_mnemonic::Symbol; verbose = :low)

    # Unfiltered population data
    population_recorded = population_data[:,[:date, population_mnemonic]]

    # Make sure first period of unfiltered population forecast is the first forecast quarter
    last_recorded_date = population_recorded[end, :date]
    if population_forecast[1, :date] <= last_recorded_date
        last_recorded_ind   = find(population_forecast[:date] .== last_recorded_date)[1]
        population_forecast = population_forecast[(last_recorded_ind+1):end, :]
    end
    @assert subtract_quarters(population_forecast[1, :date], last_recorded_date) == 1

    # population_all: full unfiltered series (including forecast)
    population_all = if isempty(population_forecast)
        population_recorded[population_mnemonic]
    else
        pop_all = vcat(population_recorded, population_forecast)
        pop_all[population_mnemonic]
    end

    # hp filter
    population_all = convert(Array{Float64}, population_all)
    filtered_population, _ = hpfilter(population_all, 1600)

    ## Setup output dictionary
    population_data_out = DataFrame()

    ## recorded series
    n_population_forecast_obs = size(population_forecast,1)

    # dates
    population_data_out[:date] = convert(Array{Date}, population_recorded[:date])

    # filtered series (levels)
    filt_pop_recorded = filtered_population[1:end-n_population_forecast_obs]
    population_data_out[:filtered_population_recorded] = filt_pop_recorded

    # filtered growth rates
    population_data_out[:dlpopulation_recorded]          = difflog(population_recorded[population_mnemonic])
    population_data_out[:dlfiltered_population_recorded] = difflog(filt_pop_recorded)

    ## forecasts
    population_forecast_out = DataFrame()
    if n_population_forecast_obs > 0
        # dates
        population_forecast_out[:date] = convert(Array{Date}, population_forecast[:date])

        # return filtered series (levels), filtered forecast growth rates, filtered data growth rates
        filt_pop_fcast = filtered_population[end-n_population_forecast_obs:end]
        population_forecast_out[:filtered_population_forecast]   = filt_pop_fcast[2:end]
        population_forecast_out[:dlpopulation_forecast]          = difflog(population_forecast[population_mnemonic])
        population_forecast_out[:dlfiltered_population_forecast] = difflog(filt_pop_fcast)[2:end]
    end

    return population_data_out, population_forecast_out
end
