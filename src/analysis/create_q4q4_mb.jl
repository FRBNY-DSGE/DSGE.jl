"""
```
create_q4q4_mb(mb)
```

Given a `MeansBands` object with `product in [:hist4q, :histforecast4q, :forecast4q, :bddforecast4q]`,
returns a `MeansBands` object with only Q4 observations.
"""
function create_q4q4_mb(mb::MeansBands)
    prod = mb.metadata[:product]
    if !(prod in [:hist4q, :histforecast4q, :forecast4q, :bddforecast4q, :bddhistforecast4q])
        error("Must pass in MeansBands with product in [:hist4q, :histforecast4q, :forecast4q, :bddforecast4q, :bddhistforecast4q] to create_forecastq4q4")
    end

    isq4(date::Date) = Dates.quarterofyear(date) == 4

    # Update metadata
    q4_metadata  = deepcopy(mb.metadata)
    q4_dates     = sort(Base.filter(isq4, collect(keys(mb.metadata[:date_inds]))))
    q4_metadata[:date_inds] = Dict(d::Date => i::Int for (i, d) in enumerate(q4_dates))
    q4_metadata[:product]   = Symbol(replace(string(prod), r"(.)4q" => s"\g<1>q4q4")) # Replace "4q" with "q4q4" in product

    # Index out Q4 observations
    q4_means = mb.means[isq4.(mb.means[!, :date]), :]

    q4_bands = Dict{Symbol, DataFrame}()
    for var in keys(mb.bands)
        dict = mb.bands[var]
        q4_bands[var] = dict[isq4.(dict[!, :date]), :]
    end

    # Return new MeansBands object
    return MeansBands(q4_metadata, q4_means, q4_bands)
end
