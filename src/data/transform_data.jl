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
(`date_forecast_start(m)` as if it were data. We also only apply transformations
for the observables given in `cond_full_names(m)` or `cond_semi_names(m)`.
"""
function transform_data(m::AbstractModel, levels::DataFrame; cond_type::Symbol = :none, verbose::Symbol = :low)

    n_obs, _ = size(levels)

    # Step 1: HP filter population forecasts, if they're being used

    population_forecast_file = if use_population_forecast(m)
        inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
    else
        ""
    end

    population_mnemonic = parse_population_mnemonic(m)[1]

    population_data = transform_population_data(levels, population_mnemonic,
                                            population_forecast_file = population_forecast_file)

    levels[:filtered_population]          = population_data[:filtered_population_recorded]
    levels[:filtered_population_growth]   = population_data[:dlfiltered_population_recorded]
    levels[:unfiltered_population_growth] = population_data[:dlpopulation_recorded]

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

    # NaN out observables not used for (semi)conditional forecasts
    if cond_type in [:semi, :full]
        cond_names = if cond_type == :semi
            cond_semi_names(m)
        elseif cond_type == :full
            cond_full_names(m)
        end

        cond_names_nan = setdiff(names(transformed), [cond_names; :date])
        T = eltype(transformed[:, cond_names_nan])
        transformed[transformed[:, :date] .>= date_forecast_start(m), cond_names_nan] = convert(T, NaN)
    end

    return transformed
end

function collect_data_transforms(m; direction=:fwd)

    data_transforms = OrderedDict{Symbol,Function}()

    # Parse vector of observable mappings into data_transforms dictionary
    for obs in keys(m.observable_mappings)
        data_transforms[obs] = getfield(m.observable_mappings[obs], symbol(string(direction) * "_transform"))
    end

    data_transforms
end

"""
```
load_population_data{S<:AbstractString}(population_data::DataFrame,
                                                 population_mnemonic::Symbol;
                                                 population_forecast_file::S = "")
```

Load, filter, and compute growth rates from population data in levels. Optionally do the same for forecasts.

### Inputs
- `population_data`: pre-loaded DataFrame containing the columns `:date` and `population_mnemonic`
- `population_mnemonic`: column name for population series in `population_data`

### Keyword Arguments
- `population_forecast_file`: filepath to a properly-formatted population forecast file.

### Output
A dictionary containing the following keys:

- `:filtered_population_recorded`: HP-filtered historical population series (levels)
- `:dlfiltered_population_recorded`: HP-filtered historical population series (growth rates)
- `:dlpopulation_recorded`: Non-filtered historical population series (growth rates)
- `:filtered_population_forecast`: HP-filtered population forecast series (levels)
- `:dlfiltered_population_forecast`: HP-filtered population forecast series (growth rates)
- `:dlpopulation_forecast`: Non-filtered historical population series (growth rates)

Note: the r"*forecast" fields will be empty if population_forecast_file is not provided.
"""
function transform_population_data{S<:AbstractString}(population_data::DataFrame,
                                                 population_mnemonic::Symbol;
                                                 population_forecast_file::S = "")

    # load historical, unfiltered population series
    population_recorded = levels[:,[:date, population_mnemonic]]

    # population_all: full unfiltered series (including forecast)
    # population_forecast: unfiltered population forecast series

    population_all, population_forecast = if !isempty(population_forecast_file)

        @assert isfile(population_forecast_file) "$population_forecast_file does not exist."

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Loading population forecast...")
        end

        # load population forecast
        pop_forecast = readtable(population_forecast_file)

        rename!(pop_forecast, :POPULATION,  population_mnemonic)
        DSGE.na2nan!(pop_forecast)
        DSGE.format_dates!(:date, pop_forecast)

        # for conditional data, start "forecast" one period later
        # (first real forecast period treated as data)
        if cond_type in [:semi, :full]
            pop_forecast = pop_forecast[2:end, :]
        end

        # use our "real" series as current value
        pop_all = vcat(population_recorded, pop_forecast[2:end, :])

        # return values
        pop_all[population_mnemonic],
        pop_forecast[population_mnemonic]
    else
        population_recorded[:,population_mnemonic], []
    end

    n_population_forecast_obs = length(population_forecast[population_mnemonic]),

    # hp filter
    population_all = convert(Array, population_all)
    filtered_population, _ = hpfilter(population_all, 1600)

    ## recorded series
    # filtered series (levels)
    filtered_population_recorded = filtered_population[1:end-n_population_forecast_obs+1]

    # filtered growth rates
    dlpopulation_recorded          = difflog(population_recorded[population_mnemonic])
    dlfiltered_population_recorded = difflog(filtered_population_recorded)

    ## forecasts
    filtered_population_forecast, dlpopulation_forecast, dlfiltered_population_forecast =
        if n_population_forecast_obs > 0
            # filtered series (levels)
            filtered_population[end-n_population_forecast_obs+1:end],

            # filtered growth rates
            difflog(population_forecast[population_mnemonic]),
            difflog(filtered_population_forecast)
        else
            [],[],[]
        end

    T = eltype(filtered_population_recorded)
    out = Dict{Symbol, Array{T}}(
        :filtered_population_recorded   => filtered_population_recorded,
        :dlfiltered_population_recorded => dlfiltered_population_recorded,
        :dlpopulation_recorded          => dlpopulation_recorded,
        :filtered_population_forecast   => filtered_population_forecast,
        :dlfiltered_population_forecast => dlfiltered_population_forecast,
        :dlpopulation_forecast          => dlpopulation_forecast
    )

    out
end