"""
```
type MeansBands
```

Stores the means and bands of results for a particular set of outputs from the forecast step.

Specifically, forecasts can be made for any element in the Cartesian product of 4 sets:

1. `input_type`: some subset of the parameter draws from the estimation
   step. See `forecast_one` for all possible options.

2. `cond_type`: conditional type. See `forecast_one` for all possible options.

3. *product*: a particular result computed in the forecast. This could be one of
   the following:

```
  - `hist`: smoothed histories
  - `forecast`: forecasted values
  - `shockdec`: shock decompositions
  - `irf`: impulse responses
```

4. variable *class*: the category in which a particular variable, like `:y_t`,
   falls. Options are:

```
  - `state`: state (from `m.endogenous_states` or `m.endogenous_states_augmented`)
  - `obs`: observable (from `m.observables`)
  - `pseudo`: pseudoobservable (from `pseudo_measurement` equation)
  - `shock`: shock (from `m.exogenous_shocks`)
```

Note that the Cartesian product (product x class) is the set of options for
`output_vars` in the `forecast_one` function signature.

### Fields

- `metadata::Dict{Symbol,Any}`: Contains metadata keeping track of the
  `input_type`, `cond_type`, product (history, forecast, shockdec,
  etc), and variable class (observable, pseudoobservable, state, etc)
  stored in this `MeansBands` structure.
- `means::DataFrame`: a `DataFrame` of the mean of the time series
- `bands::Dict{Symbol,DataFrame}`: a `Dict` mapping variable names to
  `DataFrame`s containing confidence bands for each variable. See
  `find_density_bands` for more information.
"""
type MeansBands
    metadata::Dict{Symbol,Any}
    means::DataFrame
    bands::Dict{Symbol,DataFrame}

    function MeansBands(key, means, bands)

        if !isempty(bands)

            # assert that means and bands fields have the same keys (provide info for same products)
            @assert sort(setdiff(names(means),[:date])) == sort(collect(keys(bands)))

            # check to make sure that # of periods in all dataframes are the same
            n_periods_means = size(means,1)
            for df in values(bands)
                n_periods_bands = size(df,1)
                @assert(n_periods_means == n_periods_bands,
                        "means and bands must have same number of periods")
            end
        end

        new(key, means, bands)
    end
end

