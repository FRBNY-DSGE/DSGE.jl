"""
```
reverse_transform(m, untransformed, class; verbose = :low)

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
function reverse_transform(m::AbstractModel, untransformed::DataFrame, class::Symbol;
                           verbose::Symbol = :low)
    # Dates
    @assert (:date in names(untransformed)) "untransformed must have a date column"
    date_list  = convert(Vector{Date}, untransformed[:date])
    start_date = date_list[1]
    end_date   = date_list[end]

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
    population_series = DSGE.get_population_series(:population_growth, population_data,
                                              population_forecast, start_date, end_date)

    # Calculate y0 index
    y0_date  = iterate_quarters(start_date, -1)
    y0_index = findfirst(data[:date], y0_date)

    # Apply reverse transform
    transformed = DataFrame()
    transformed[:date] = untransformed[:date]

    var_names = setdiff(names(untransformed), [:date])
    for var in var_names
        rev_transform = dict[var].rev_transform
        y = convert(Vector{Float64}, untransformed[var])
        y0 = if y0_index > 0 && class == :obs
            data[y0_index, var]
        else
            NaN
        end

        transformed[var] = reverse_transform(y, rev_transform;
                               y0 = y0, pop_growth = population_series)
    end
    return transformed
end

function reverse_transform{T<:AbstractFloat}(y, rev_transform::Function;
                                             y0::T = NaN,
                                             pop_growth::Vector{T} = Vector{T}(),
                                             y0_index::Int = -1)

    if rev_transform in [loggrowthtopct_annualized_percapita]
        rev_transform(y, pop_growth)
    elseif rev_transform in [logleveltopct_annualized_percapita]
        rev_transform(y, y0, pop_growth)
    else
        rev_transform(y)
    end
end
