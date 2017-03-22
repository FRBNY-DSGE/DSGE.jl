"""
```
reverse_transform(m, untransformed, class, product; verbose = :low)

reverse_transform(series, rev_transform; data_series = [],
    population_series = [], y0_index = -1))
```

The first method takes a dated `DataFrame` of series in model units. It uses the
model object to load data and population data, as well as look up the
appropriate reverse transformations for each series. It returns a `DataFrame` of
the transformed series.

The second method takes a single series and reverse transformation, applies it,
and returns a vector.
"""
function reverse_transform(m::AbstractModel, untransformed::DataFrame, class::Symbol, product::Symbol;
                           verbose::Symbol = :low)
    # Dates
    @assert (:date in names(untransformed)) "untransformed must have a date column"
    date_list = convert(Vector{Date}, untransformed[:date])

    # Get mapping from variable name to Observable or PseudoObservable instance
    if class == :obs
        dict = m.observable_mappings
    elseif class == :pseudo
        dict, _ = pseudo_measurement(m)
    else
        error("Class $class does not have reverse transformations")
    end

    # Load data
    data = load_data(m; verbose = verbose)

    # Prepare population series
    population_mnemonic = DSGE.parse_population_mnemonic(m)[1]
    vint = data_vintage(m)
    population_data_file     = inpath(m, "data", "population_data_levels_$vint.csv")
    population_forecast_file = inpath(m, "data", "population_forecast_$vint.csv")
    population_data, population_forecast =
        DSGE.load_population_growth(population_data_file, population_forecast_file,
                               get(population_mnemonic); verbose = verbose)
    population_series = DSGE.get_mb_population_series(product, :population_growth, population_data,
                                                 population_forecast, date_list)

    # Apply reverse transform
    transformed = DataFrame()
    transformed[:date] = untransformed[:date]
    y0_index = DSGE.get_y0_index(m, product)

    var_names = setdiff(names(untransformed), [:date])
    for var in var_names
        series = convert(Vector{Float64}, untransformed[var])
        rev_transform = dict[var].rev_transform
        data_series = convert(Vector{Float64}, data[var])

        transformed[var] = reverse_transform(series, rev_transform;
                               data_series = data_series,
                               population_series = population_series,
                               y0_index = y0_index)
    end
    return transformed
end

function reverse_transform{T<:AbstractFloat}(series, rev_transform::Function;
                                             data_series::Vector{T} = Vector{T}(),
                                             population_series::Vector{T} = Vector{T}(),
                                             y0_index::Int = -1)

    if rev_transform in [logtopct_annualized_percapita]
        rev_transform(series, population_series)
    elseif rev_transform in [logleveltopct_annualized_percapita]
        rev_transform(series, data_series[y0_index], population_series)
    else
        rev_transform(series)
    end
end
