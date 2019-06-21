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

    # Step 1: HP filter (including population forecasts, if they're being used)
    population_mnemonic = parse_population_mnemonic(m)[1]
    if !isnull(population_mnemonic)
        population_forecast_levels = if use_population_forecast(m)
            read_population_forecast(m; verbose = verbose)
        else
            DataFrame()
        end

        population_data, _ = transform_population_data(levels, population_forecast_levels,
                                                       get(population_mnemonic);
                                                       verbose = verbose,
                                                       use_hpfilter = hpfilter_population(m))

        levels = join(levels, population_data, on = :date, kind = :left)
        name_maps = [s => t for (s,t) = zip([:filtered_population_recorded, :dlfiltered_population_recorded, :dlpopulation_recorded],
                [:filtered_population, :filtered_population_growth, :unfiltered_population_growth])]
        rename!(levels, name_maps)
    end

    # Step 2: apply transformations to each series
    transformed = DataFrame()
    transformed[:date] = levels[:date]

    data_transforms = collect_data_transforms(m)

    for series in keys(data_transforms)
        println(verbose, :high, "Transforming series $series...")
        f = data_transforms[series]
        transformed[series] = f(levels)
    end

    sort!(transformed, [:date])

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
- `use_hpfilter`: whether to HP filter population data and forecast. See `Output` below.
- `pad_forecast_start::Bool`: Whether you want to re-size
the population_forecast such that the first index is one quarter ahead of the last index
of population_data. Only set to false if you have manually constructed population_forecast
to artificially start a quarter earlier, so as to avoid having an unnecessary missing first entry.

### Output

Two dictionaries containing the following keys:

- `population_data_out`:
  + `:filtered_population_recorded`: HP-filtered historical population series (levels)
  + `:dlfiltered_population_recorded`: HP-filtered historical population series (growth rates)
  + `:dlpopulation_recorded`: Non-filtered historical population series (growth rates)

- `population_forecast_out`:
  + `:filtered_population_forecast`: HP-filtered population forecast series (levels)
  + `:dlfiltered_population_forecast`: HP-filtered population forecast series (growth rates)
  + `:dlpopulation_forecast`: Non-filtered historical population series (growth rates)

If `population_forecast_file` is not provided, the r\"*forecast\" fields will be
empty. If `use_hpfilter = false`, then the r\"*filtered*\" fields will be
empty.
"""
function transform_population_data(population_data::DataFrame, population_forecast::DataFrame,
                                   population_mnemonic::Symbol; verbose = :low,
                                   use_hpfilter::Bool = true,
                                   pad_forecast_start::Bool = false)

    # Unfiltered population data
    population_recorded = population_data[[:date, population_mnemonic]]

    # Make sure first period of unfiltered population forecast is the first forecast quarter
    if !isempty(population_forecast) && !pad_forecast_start
        last_recorded_date = population_recorded[end, :date]
        if population_forecast[1, :date] <= last_recorded_date
            last_recorded_ind   = findall(population_forecast[:date] .== last_recorded_date)[1]
            population_forecast = population_forecast[(last_recorded_ind+1):end, :]
        end
        @assert subtract_quarters(population_forecast[1, :date], last_recorded_date) == 1
    end

    # population_all: full unfiltered series (including forecast)
    population_all = if isempty(population_forecast)
        population_recorded[population_mnemonic]
    else
        pop_all = vcat(population_recorded, population_forecast)
        pop_all[population_mnemonic]
    end

    # HP filter
    if use_hpfilter
        population_all = convert(Array{Union{Float64, Missing}}, population_all)
        filtered_population, _ = hpfilter(population_all, 1600)
    end

    # Output dictionary for population data
    population_data_out = DataFrame()
    population_data_out[:date] = convert(Array{Date}, population_recorded[:date])
    population_data_out[:dlpopulation_recorded] = difflog(population_recorded[population_mnemonic])

    n_population_forecast_obs = size(population_forecast,1)

    if use_hpfilter
        filt_pop_recorded = filtered_population[1:end-n_population_forecast_obs]
        population_data_out[:filtered_population_recorded] = filt_pop_recorded
        population_data_out[:dlfiltered_population_recorded] = difflog(filt_pop_recorded)
    end

    # Output dictionary for population forecast
    population_forecast_out = DataFrame()
    if n_population_forecast_obs > 0
        population_forecast_out[:date] = convert(Array{Date}, population_forecast[:date])
        population_forecast_out[:dlpopulation_forecast] = difflog(population_forecast[population_mnemonic])

        if use_hpfilter
            filt_pop_fcast = filtered_population[end-n_population_forecast_obs:end]
            population_forecast_out[:filtered_population_forecast]   = filt_pop_fcast[2:end]
            population_forecast_out[:dlfiltered_population_forecast] = difflog(filt_pop_fcast)[2:end]
        end
    end

    return population_data_out, population_forecast_out
end
