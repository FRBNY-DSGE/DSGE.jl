"""
```
reverse_transform(m, untransformed, class; fourquarter = false, verbose = :low)

reverse_transform(y, rev_transform; fourquarter = false, y0 = NaN, y0s = [],
    pop_growth = [])
```

The first method takes a dated `DataFrame` of series in model units. It uses the
model object to load population data, as well as look up the appropriate reverse
transformations for each series. It returns a `DataFrame` of the transformed
series.

The second method takes a single series and reverse transformation, applies it,
and returns a vector.

## Inputs
- `m::AbstractModel`
- `untransformed`: `nperiods` x `nseries` `Matrix` or `DataFrame` of series in model units.

## Keyword arguments
- `fourquarter::Bool`: produce four-quarter
   results instead of annualized quarter-to-quarter results
"""
function reverse_transform(m::AbstractModel, untransformed::Matrix, start_date::Date,
                           var_names::Vector{Symbol}, class::Symbol;
                           fourquarter::Bool = false,
                           verbose::Symbol = :low)
    df = DataFrame()
    nperiods = size(untransformed, 2)
    end_date = iterate_quarters(start_date, nperiods - 1)
    df[:date] = quarter_range(start_date, end_date)

    for (ind, var) in enumerate(var_names)
        series = squeeze(untransformed[ind, :], 1)
        df[var] = series
    end

    reverse_transform(m, df, class; fourquarter = fourquarter, verbose = verbose)
end

function reverse_transform(m::AbstractModel, untransformed::DataFrame, class::Symbol;
                           fourquarter::Bool = false, verbose::Symbol = :low)
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

    # Prepare population series
    population_mnemonic = parse_population_mnemonic(m)[1]
    vint = data_vintage(m)
    population_data_file     = inpath(m, "data", "population_data_levels_$vint.csv")
    population_forecast_file = inpath(m, "data", "population_forecast_$vint.csv")
    population_data, population_forecast =
        load_population_growth(population_data_file, population_forecast_file,
                               get(population_mnemonic); verbose = verbose)
    population_series = get_population_series(:population_growth, population_data,
                                              population_forecast, start_date, end_date)

    # Apply reverse transform
    transformed = DataFrame()
    transformed[:date] = untransformed[:date]

    var_names = setdiff(names(untransformed), [:date])
    for var in var_names
        y = convert(Vector{Float64}, untransformed[var])
        rev_transform = dict[var].rev_transform
        if fourquarter
            rev_transform = get_transform4q(rev_transform)
        end
        transformed[var] = reverse_transform(y, rev_transform; fourquarter = fourquarter,
                                             pop_growth = population_series)
    end
    return transformed
end

function reverse_transform{T<:AbstractFloat}(y::Array{T}, rev_transform::Function;
                                             fourquarter::Bool = false,
                                             y0::T = NaN, y0s::Vector{T} = Vector{T}(),
                                             pop_growth::Vector{T} = Vector{T}())
    if fourquarter
        if rev_transform in [loggrowthtopct_4q_percapita, loggrowthtopct_4q]
            # Sum growth rates y_{t-3}, y_{t-2}, y_{t-1}, and y_t
            y0s = isempty(y0s) ? fill(NaN, 3) : y0s
            if rev_transform == loggrowthtopct_4q_percapita
                rev_transform(y, y0s, pop_growth)
            else
                rev_transform(y, y0s)
            end
        elseif rev_transform in [logleveltopct_4q_percapita, logleveltopct_4q]
            # Divide log levels y_t by y_{t-4}
            y0s = isempty(y0s) ? fill(NaN, 4) : y0s
            if rev_transform == logleveltopct_4q_percapita
                rev_transform(y, y0s, pop_growth)
            else
                rev_transform(y, y0s)
            end
        elseif rev_transform in [quartertoannual, identity]
            rev_transform(y)
        else
            error("Invalid 4-quarter reverse transform: $rev_transform")
        end
    else
        if rev_transform in [loggrowthtopct_annualized_percapita]
            rev_transform(y, pop_growth)
        elseif rev_transform in [logleveltopct_annualized_percapita]
            rev_transform(y, y0, pop_growth)
        elseif rev_transform in [loggrowthtopct_annualized, logleveltopct_annualized, quartertoannual, identity]
            rev_transform(y)
        else
            error("Invalid reverse transform: $rev_transform")
        end
    end
end
