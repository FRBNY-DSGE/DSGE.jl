"""
```
mutable struct MeansBands
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
  - `shockdecseq`: shock decompositions
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
mutable struct MeansBands
    metadata::Dict{Symbol,Any}
    means::DataFrame
    bands::Dict{Symbol,DataFrame}

    function MeansBands(key, means, bands)

        if !isempty(bands)

            # assert that means and bands fields have the same keys (provide info for same products)
            @assert sort(setdiff(propertynames(means),[:date])) == sort(collect(keys(bands)))

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

# A dummy MeansBands object
function MeansBands()
    metadata   = Dict(:class => :none, :product => :none,
                      :cond_type => :none, :para => :none,
                      :date_inds => OrderedDict{Date, Int}(Date(1) => 1),
                      :indices => Dict{Symbol, Int}(:none => 1))

    means = DataFrame(date = [Dates.Date(0)], none = [0.0])
    bands = Dict{Symbol,DataFrame}(:none => DataFrame(date = [Dates.Date(0)]))

    MeansBands(metadata, means, bands)
end

function Base.show(io::IO, mb::MeansBands)
    @printf io "MeansBands\n"
    @printf io "  class: %s\n"   get_class(mb)
    @printf io "  product: %s\n" get_product(mb)
    @printf io "  cond: %s\n"    get_cond_type(mb)
    @printf io "  para: %s\n"    get_para(mb)
    if mb.metadata[:product] != :trend && mb.metadata[:product] != :irf
        if isempty(mb.metadata[:date_inds])
            @printf io "  dates: []\n"
        else
            @printf io "  dates: %s - %s\n" startdate_means(mb) enddate_means(mb)
        end
    end
    @printf io "  # of variables: %s\n" n_vars_means(mb)
    @printf io "  bands: %s\n" which_density_bands(mb, uniquify=true)
    if haskey(mb.metadata, :scenario_key)
        @printf io "  scenario: %s\n" get_scenario_key(mb)
    end
end

"""
```
Base.isempty(mb::MeansBands)
```

Returns whether the `mb` object in question is a dummy.
"""
function Base.isempty(mb::MeansBands)

    return get_class(mb) == :none && get_product(mb) == :none &&
        startdate_means(mb) == Dates.Date(0) &&
        collect(propertynames(mb.means)) ==  [:date, :none] &&
        collect(keys(mb.bands)) ==  [:none]
end

"""
```
Base.cat(mb1::MeansBands, mb2::MeansBands; out_product = Symbol(),
    forecast_string = "")
```

Concatenate 2 compatible `MeansBands` objects together by date. 2
`MeansBands` objects are defined to be compatible if the class of
variables is the same, the conditional type is the same, and the
input_type for the forecast used to create the two `MeansBands` object
is the same. Furthermore, we require that the dates covered by each
`MeansBands` object form a continguous interval.

### Inputs

- `mb1::MeansBands`
- `mb2::MeansBands`

Note that the dates in `mb1` should come chronologically first, with the dates
in `mb2` beginning 1 period after the final period in `mb1`.

### Keyword Arguments

- `out_product::Symbol`: desired product of the resulting concatenated
`MeansBands` object. This argument is required if `mb1` and `mb2` do not
represent a history and a forecast respectively
- `forecast_string::String`: desired `forecast_string` of the resulting
  concatenated `MeansBands`. This argument is recommended (but not required) if
  the `forecast_string`s of `mb1` and `mb2` do not match. If this is the case
  but `forecast_string` is not provided, `mb1`'s `forecast_string` will be used.
"""
function Base.cat(mb1::MeansBands, mb2::MeansBands;
                  out_product::Symbol = Symbol(),
                  forecast_string::String = "")

    # If either mb1 or mb2 is empty, return just the other one
    if isempty(mb1)
        return mb2
    elseif isempty(mb2)
        return mb1
    end

    # Assert class, cond type, para, and scenario info are the same
    @assert get_class(mb1) == get_class(mb2)
    @assert get_cond_type(mb1) == get_cond_type(mb2)
    @assert get_para(mb1) == get_para(mb2)
    @assert haskey(mb1.metadata, :scenario_key) == haskey(mb2.metadata, :scenario_key) ==
        haskey(mb1.metadata, :scenario_vint) == haskey(mb2.metadata, :scenario_vint)
    if haskey(mb1.metadata, :scenario_key)
        @assert mb1.metadata[:scenario_key] == mb2.metadata[:scenario_key]
        @assert mb1.metadata[:scenario_vint] == mb2.metadata[:scenario_vint]
    end

    # Assert dates are contiguous
    last_mb1_date  = enddate_means(mb1)
    first_mb2_date = startdate_means(mb2)
#    @assert iterate_quarters(last_mb1_date, 1) == first_mb2_date

    # compute means field
    means = vcat(mb1.means, mb2.means)
    # na2nan!(means) # Removed b/c DataFrames uses missing instead of NA, and missings are already handled

    # compute bands field
    bands = Dict{Symbol, DataFrame}()

    mb1vars = collect(keys(mb1.bands))
    mb2vars = collect(keys(mb2.bands))
    nperiods_mb1 = length(mb1.metadata[:date_inds])
    nperiods_mb2 = length(mb2.metadata[:date_inds])

    bothvars = intersect(mb1vars, mb2vars)
    for var in setdiff(union(keys(mb1.bands), keys(mb2.bands)), [:date])
        bands[var] = if var in bothvars
            vcat(mb1.bands[var], mb2.bands[var])
        elseif var in setdiff(mb1vars, mb2vars)
            vcat(mb1vars[var], fill(NaN, nperiods_mb2))
        else
            vcat(fill(NaN, nperiods_mb1), mb2vars[var])
        end
        # na2nan!(bands[var]) # Removed b/c DataFrames uses missing instead of NA, and missings are already handled
    end

    # compute metadata
    # product
    mb1_product = get_product(mb1)
    mb2_product = get_product(mb2)
    product = if mb1_product == :hist
        if mb2_product == :forecast
            :histforecast
        elseif mb2_product == :bddforecast
            :bddforecast
        end
    elseif mb1_product == :histut
        if mb2_product == :forecastut
            :histforecastut
        elseif mb2_product == :bddforecast
            :bddforecastut
        end
    elseif mb1_product == :hist4q
        if mb2_product == :forecast4q
            :histforecast4q
        elseif mb2_product == :bddforecast4q
            :bddhistforecast4q
        end
    elseif mb1_product == mb2_product
        mb1_product
    else
        @assert !isempty(out_product) "Please supply a product name for the output MeansBands"
        out_product
    end

    # date indices
    date_indices = Dict(d::Dates.Date => i::Int for (i, d) in enumerate(means[!, :date]))

    # variable indices
    indices = Dict(var::Symbol => i::Int for (i, var) in enumerate(propertynames(means)))

    # forecast string
    if isempty(forecast_string) && (mb1.metadata[:forecast_string] != mb2.metadata[:forecast_string])
        @warn "No forecast_string provided: using $(mb1.metadata[:forecast_string])"
    end
    forecast_string = mb1.metadata[:forecast_string]

    metadata = Dict{Symbol, Any}(
                   :para            => get_para(mb1),
                   :cond_type       => get_cond_type(mb1),
                   :class           => get_class(mb1),
                   :product         => product,
                   :indices         => indices,
                   :forecast_string => forecast_string,
                   :date_inds       => sort(date_indices, by = x -> date_indices[x]))
    if haskey(mb1.metadata, :scenario_key)
        metadata[:scenario_key] = mb1.metadata[:scenario_key]
        metadata[:scenario_vint] = mb1.metadata[:scenario_vint]
    end

    # construct the new MeansBands object and return
    MeansBands(metadata, means, bands)
end


###################################
## METADATA
###################################

"""
```
get_class(mb::MeansBands)
```
Returns the class of variables contained this `MeansBands` object (observables, pseudoobservables)
"""
get_class(mb::MeansBands) = mb.metadata[:class]

"""
```
get_product(mb::MeansBands)
```
Returns the product stored in this `MeansBands` object (history, forecast, shock decompostion, etc).
"""
get_product(mb::MeansBands) = mb.metadata[:product]

"""
```
get_cond_type(mb::MeansBands)
```
Returns the conditional type of the forecast used to create this `MeansBands` object.
"""
get_cond_type(mb::MeansBands) = mb.metadata[:cond_type]

"""
```
get_para(mb::MeansBands)
```
Returns the `input_type` from the forecast used to create this `MeansBands` object.
"""
get_para(mb::MeansBands) = mb.metadata[:para]

"""
```
get_shocks(mb::MeansBands)
```

Returns a list of shock names that are used for the shock
decomposition stored in a shock decomposition or irf MeansBands object `mb`.
"""
function get_shocks(mb::MeansBands)
    @assert get_product(mb) in [:shockdec, :irf, :shockdecseq] "Function only for shockdec, shockdecseq, or irf MeansBands objects"
    varshocks = setdiff(propertynames(mb.means), [:date])
    unique(map(x -> Symbol(split(string(x), DSGE_SHOCKDEC_DELIM)[2]), varshocks))
end

"""
```
parse_mb_colname(s::Symbol)
```

`MeansBands` column names are saved in the format
`\$var\$DSGE_SHOCKDEC_DELIM\$shock`. `parse_mb_colname` returns (`var`,
`shock`).
"""
function parse_mb_colname(s::Symbol)
    map(Symbol, split(string(s), DSGE_SHOCKDEC_DELIM))
end

"""
```
get_variables(mb::MeansBands)
```

Returns a list of variable names that are used for the shock
decomposition stored in a shock decomposition or irf MeansBands object `mb`.
"""
function get_variables(mb::MeansBands)
    @assert get_product(mb) in [:shockdec, :irf, :shockdecseq] "Function only for shockdec or irf MeansBands objects"
    varshocks = setdiff(propertynames(mb.means), [:date])
    unique(map(x -> Symbol(split(string(x), DSGE_SHOCKDEC_DELIM)[1]), varshocks))
end

"""
```
get_scenario_key(mb::MeansBands)
```

If `mb` is an alternative scenario `MeansBands`, return the scenario key.
"""
function get_scenario_key(mb::MeansBands)
    @assert haskey(mb.metadata, :scenario_key) "Function only for scenario MeansBands objects"
    mb.metadata[:scenario_key]
end


###################################
## MEANS
###################################

"""
```
n_vars_means(mb::MeansBands)
````

Get number of variables (`:y_t`, `:OutputGap`, etc) in `mb.means`.
"""
function n_vars_means(mb::MeansBands)
    length(get_vars_means(mb))
end

"""
```
get_vars_means(mb::MeansBands)
````

Get variables (`:y_t`, `:OutputGap`, etc) in `mb.means`. Note that
`mb.metadata[:indices]` is an `OrderedDict`, so the keys will be in the correct
order.
"""
function get_vars_means(mb::MeansBands)
    collect(keys(mb.metadata[:indices]))
end


"""
```
n_periods_means(mb::MeansBands)
```

Get number of periods in `mb.means`.
"""
n_periods_means(mb::MeansBands) = size(mb.means,1)

"""
```
startdate_means(mb::MeansBands)
```

Get first period in`mb.means`. Assumes `mb.means[product]` is already sorted by
date.
"""
startdate_means(mb::MeansBands) = mb.means[!, :date][1]

"""
```
enddate_means(mb::MeansBands)
```

Get last period for which `mb` stores means. Assumes `mb.means[product]` is
already sorted by date.
"""
enddate_means(mb::MeansBands) = mb.means[!, :date][end]


"""
```
get_shockdec_means(mb::MeansBands, var::Symbol;
    shocks::Vector{Symbol} = Vector{Symbol}())
```

Return the mean value of each shock requested in the shock decomposition of a
particular variable. If `shocks` is empty, returns all shocks.
"""
function get_shockdec_means(mb::MeansBands, var::Symbol; shocks::Vector{Symbol} = Vector{Symbol}())

    # Extract the subset of columns relating to the variable `var` and the shocks listed in `shocks.`
    # If `shocks` not provided, give all the shocks
    var_cols = collect(propertynames(mb.means))[findall([contains(string(col), string(var)) for col in propertynames(mb.means)])]
    if !isempty(shocks)
        var_cols = [col -> contains(string(col), string(shock)) ? col : nothing for shock in shocks]
    end

    # Make a new DataFrame with column the column names
    out = DataFrame()
    for col in var_cols
        shockname = split(string(col), DSGE_SHOCKDEC_DELIM)[2]
        out[Symbol(shockname)] = mb.means[col]
    end

    return out
end



###################################
## BANDS
###################################

"""
```
n_vars_bands(mb::MeansBands)
```

Get number of variables (`:y_t`, `:OutputGap`, etc) for which `mb`
stores bands for the specified `product` (`hist`, `forecast`, `shockdec`, etc).
"""
n_vars_bands(mb::MeansBands) = length(mb.bands)


"""
```
n_periods_bands(mb::MeansBands)
```

Get number of periods for which `mb` stores bands for the specified
product` (`hist`, `forecast`, `shockdec`, etc).
"""
function n_periods_bands(mb::MeansBands)
    size(mb.bands[collect(keys(mb.bands))[1]],1)
end

"""
```
startdate_bands(mb::MeansBands)
```

Get first period for which `mb` stores bands. Assumes `mb.bands` is already sorted by date.
Also assumes the startdate is same for all observables in MeansBands
"""
startdate_bands(mb::MeansBands) = if isempty(mb)
    Date(0000, 1, 1)
else
    first(mb.bands)[2][:, :date][1]
end

"""
```
enddate_bands(mb::MeansBands)
```

Get last period in `mb.bands`. Assumes `mb.bands` is already sorted by date.
Also assumes the startdate is same for all observables in MeansBands
"""
enddate_bands(mb::MeansBands) = if isempty(mb)
    Date(0000, 1, 1)
else
    first(mb.bands)[2][:, :date][end]
end

"""
```
which_density_bands(mb, uniquify=false)
```

Return a list of the bands stored in mb.bands. If `uniquify = true`,
strips \"upper\" and \"lower\" band tags and returns unique list of percentage values.
"""
function which_density_bands(mb::MeansBands; uniquify=false, ordered=true)

    if isempty(mb.bands)
        return String[]
    else
        # extract one of the keys in mb.bands
        var  = collect(keys(mb.bands))[1]

        # get all the columns in the corresponding dataframe that aren't dates
        strs = map(string,propertynames(mb.bands[var]))
        strs = setdiff(strs, ["date"])

        lowers = strs[map(occursin, repeat([r"LB"], outer=length(strs)), strs)]
        uppers = strs[map(occursin, repeat([r"UB"], outer=length(strs)), strs)]

        # sort
        if ordered
            sort!(lowers, rev=true)
            sort!(uppers)
        end

        # return both upper and lower bands, or just percents, as desired
        strs = if uniquify
            sort([convert(String, split(x, " ")[1]) for x in lowers])
        else
            [lowers; uppers]
        end
        return strs
    end
end


"""
```
get_shockdec_bands(mb, var; shocks = Vector{Symbol}(), bands = Vector{Symbol}())
```

Return a `Dict{Symbol,DataFrame}` mapping shock names to bands for a particular
variable.

### Inputs

- `mb::MeansBands`
- `var::Symbol`: the variable of interest (eg the state `:y_t`, or observable
  `:obs_hours`)

### Keyword Arguments

- `shocks::Vector{Symbol}`: subset of shock names for which to return bands. If
  empty, `get_shockdec_bands` returns all bands
- `bands::Vector{Symbol}`: subset of bands stored in the DataFrames of
  `mb.bands` to return
"""
function get_shockdec_bands(mb::MeansBands, var::Symbol;
                            shocks::Vector{Symbol} = Vector{Symbol}(),
                            bands::Vector{Symbol} = Vector{Symbol}())

    @assert get_product(mb) == :shockdec || get_product(mb) == :shockdecseq

    # Extract the subset of columns relating to the variable `var` and the shocks listed in `shocks.`
    # If `shocks` not provided, give all the shocks
    var_cols = collect(keys(mb.bands))[findall([contains(string(col), string(var)) for col in keys(mb.bands)])]
    if !isempty(shocks)
        var_cols = [col -> contains(string(col), string(shock)) ? col : nothing for shock in shocks]
    end

    # Extract the subset of bands we want to return. Return all bands if `bands` not provided.
    bands_keys = if isempty(bands)
        propertynames(mb.bands[var_cols[1]])
    else
        [[Symbol("$(100x)% LB") for x in bands]; [Symbol("$(100x)% UB") for x in bands]]
    end

    # Make a new dictionary mapping shock names to bands
    out = Dict{Symbol, DataFrame}()
    for col in var_cols
        shockname = parse_mb_colname(col)[2]
        out[shockname] = mb.bands[col][bands_keys]
    end

    return out
end


################################################
## EXTRACTING VARIABLES
################################################

"""
```
prepare_meansbands_table_timeseries(mb, var)
```

Returns a `DataFrame` of means and bands for a particular time series variable
(either `hist` or `forecast` of some type). Columns are sorted such that the
bands are ordered from smallest to largest, and the means are on the far
right. For example, a `MeansBands` containing 50\\% and 68\\% bands would be
ordered as follows: [68\\% lower, 50\\% lower, 50\\% upper, 68\\% upper, mean].

### Inputs

- `mb::MeansBands`: time-series MeansBands object
- `var::Symbol`: an economic variable stored in `mb`. If `mb` stores
  observables, `var` would be an element of `propertynames(m.observables)`. If
  it stores pseudo-observables, `var` would be the name of a
  pseudo-observable defined in the pseudo-measurement equation.

### Keyword Arguments

- `bands_pcts::Vector{String}`: vector of (uniquified) band percentiles to
  include in the table
"""
function prepare_meansbands_table_timeseries(mb::MeansBands, var::Symbol;
                                             bands_pcts::Vector{String} = which_density_bands(mb, uniquify = true))

    @assert get_product(mb) in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                                :bddforecast, :bddforecastut, :bddforecast4q,
                                :histforecast, :histforecastut, :histforecast4q,
                                :bddhistforecast, :bddhistforecastut, :bddhistforecast4q,
                                :trend, :dettrend] "prepare_meansbands_table_timeseries can only be used for time-series products"
    @assert var in get_vars_means(mb) "$var is not stored in this MeansBands object"

    # Get bands
    uppers = sort!([pct * " UB" for pct in bands_pcts], rev = true)
    lowers = sort!([pct * " LB" for pct in bands_pcts])
    my_bands = map(Symbol, vcat(lowers, uppers))

    # Extract this variable from Means and bands
    means = mb.means[:, [:date, var]]
    bands = mb.bands[var][:, [:date; my_bands]]

    # Join so mean is on far right and date is on far left
    df = isdefined(DataFrames, :innerjoin) ? innerjoin(bands, means, on = :date) : join(bands, means, on = :date)
    rename!(df, var => Symbol("mean"))

    return df
end

"""
```
prepare_means_table_shockdec(mb_shockdec, mb_trend, mb_dettrend, var;
    shocks = get_shocks(mb_shockdec), mb_forecast = MeansBands(),
    mb_hist = MeansBands(), detexify = true, groups = [],
    states_trend_shock = false)
```

Returns a `DataFrame` representing a detrended shock decompostion for
the variable `var`. The columns of this dataframe represent the
contributions of each shock in `shocks` (or all shocks, if the keyword
argument is omitted) and the deterministic trend.

### Inputs

- `mb_shockdec::MeansBands`: a `MeansBands` object for a shock decomposition
- `mb_trend::MeansBands`: a `MeansBands` object for a trend  product.
- `mb_dettrend::MeansBands`: a `MeansBands` object for a deterministic trend
  product.
- `var::Symbol`: name of economic variable for which to return the means and bands table

### Keyword Arguments

- `shocks::Vector{Symbol}`: If `mb` is a shock decomposition, this is
  an optional list of shocks to print to the table. If omitted, all
  shocks will be printed.
- `mb_forecast::MeansBands`: a `MeansBands` object for a forecast.
- `mb_hist::MeansBands`: a `MeansBands` object for smoothed states.
- `detexify::Bool`: whether to remove Unicode characters from shock names
- `groups::Vector{ShockGroup}`: if provided, shocks will be grouped accordingly
- `states_trend_shock::Bool`: whether to treat nonzero trends in the states
    from nonzero `CCC` as a shock rather than a trend.
"""
function prepare_means_table_shockdec(mb_shockdec::MeansBands, mb_trend::MeansBands,
                                      mb_dettrend::MeansBands, var::Symbol;
                                      shocks::Vector{Symbol} = get_shocks(mb_shockdec),
                                      mb_forecast::MeansBands = MeansBands(),
                                      mb_hist::MeansBands = MeansBands(),
                                      detexify_shocks::Bool = true,
                                      groups::Vector{ShockGroup} = ShockGroup[],
                                      trend_nostates::DataFrame = DataFrame(), df_enddate = Date(2100,12,31))

    @assert get_product(mb_shockdec) in [:shockdec, :shockdecseq] "The first argument must be a MeansBands object for a shockdec"
    @assert get_product(mb_trend)    == :trend    "The second argument must be a MeansBands object for a trend"
    @assert get_product(mb_dettrend) == :dettrend "The third argument must be a MeansBands object for a deterministic trend"

    # Get the variable-shock combinations we want to print
    varshocks = [Symbol("$var" * DSGE_SHOCKDEC_DELIM * "$shock") for shock in shocks]

    # Fetch the columns corresponding to varshocks
    df_shockdec = mb_shockdec.means[!,union([:date], varshocks)]
    df_trend    = mb_trend.means[!,[:date, var]]
    df_dettrend = mb_dettrend.means[!,[:date, var]]

    # Line up dates between trend, dettrend and shockdec
    has_ij = isdefined(DataFrames, :innerjoin)
    if has_ij
        df_shockdec = innerjoin(df_shockdec, df_trend, on = :date)
        rename!(df_shockdec, var => :trend)
        df_shockdec = innerjoin(df_shockdec, df_dettrend, on = :date)
        rename!(df_shockdec, var => :dettrend)
    else
        df_shockdec = join(df_shockdec, df_trend, on = :date, kind = :inner)
        rename!(df_shockdec, var => :trend)
        df_shockdec = join(df_shockdec, df_dettrend, on = :date, kind = :inner)
        rename!(df_shockdec, var => :dettrend)
    end

    # Add each shock's contribution and deterministic trend to output DataFrame
    df = DataFrame(date = df_shockdec[!,:date])
    for col in setdiff(propertynames(df_shockdec), [:date, :trend])
        df[!,col] = df_shockdec[!,col]
    end
    df_shockdec[!,:dettrend] = df_shockdec[!,:dettrend]

    # Rename columns to just the shock names
    map(x -> rename!(df, x => parse_mb_colname(x)[2]), setdiff(propertynames(df), [:date, :trend, :dettrend]))

    # If mb_forecast and mb_hist are passed in, add the detrended time series
    # mean of var to the table
    if !(isempty(mb_forecast) && isempty(mb_hist))
        mb_timeseries = cat(mb_hist, mb_forecast)

        # Truncate to just the dates we want
        startdate = df[1, :date]
        enddate   = min(df[end, :date], df_enddate)
        df_mean   = mb_timeseries.means[startdate .<= mb_timeseries.means[!, :date] .<= enddate, [:date, var]]

        df_shockdec = has_ij ? innerjoin(df_shockdec, df_mean, on = :date) :
            join(df_shockdec, df_mean, on = :date, kind = :inner)
        df = df[startdate .<= df[!, :date] .<= enddate, :]

        if isempty(trend_nostates)
            df[!, :detrendedMean] = df_shockdec[!, var] - df_shockdec[!, :trend]
        else
            var_trend_nostates   = trend_nostates[startdate .<= trend_nostates[!, :date] .<= enddate, var]
            if size(var_trend_nostates, 1) != size(df_shockdec, 1)
                error("The number of rows in kwarg `trend_nostates` does not match the number in df_shockdec. Check that " *
                      "the Setting :date_forecast_end matches the date used for the calculation of the shockdecs " *
                      "when constructing the `trend_nostates` DataFrame.")
            end

            var_trend_states      = df_shockdec[!, :trend] - var_trend_nostates
            df[!, :detrendedMean] = df_shockdec[!,var] - var_trend_nostates
            df[!, :StatesTrend]   = var_trend_states
        end
    end

    # Group shocks if desired
    nperiods = size(df, 1)
    v0 = zeros(nperiods)

    rm_shocks = []
    for group in groups
        # Sum shock values for each group
        shock_vectors = [df[!,shock] for shock in group.shocks]
        shock_sum = reduce(+, shock_vectors; init = v0)
        df[!,Symbol(group.name)] = shock_sum

        # Delete original (ungrouped) shocks from df
        rm_shocks = vcat(rm_shocks, setdiff(group.shocks, [Symbol(group.name)]))
        # select!(df, Not(setdiff(group.shocks, [Symbol(group.name)])))
    end
    select!(df, Not(unique(rm_shocks)))

    # Remove Unicode characters from shock names
    if detexify_shocks
        for x in setdiff(propertynames(df), [:date, :trend, :dettrend])
            x_detexed = detexify(x)
            if x != x_detexed
                rename!(df, x => x_detexed)
            end
        end
    end

    return df
end

# Helper function for plotting the trend in states separately. This function
# calculates trends without any trends in states. For example,
# rather than return the trend in states, this function will return a DataFrame of
# all zeros since the trend in states without the trend in states is just detrending the trend!
# This function is tested in test/forecast/shock_decompositions.jl
function prepare_means_table_trend_nostates(m::AbstractDSGEModel{S}, cond_type::Symbol, class::Symbol,
                                            start_date::Date, end_fcast_date::Date; annualize::Bool = true,
                                            apply_altpolicy::Bool = false) where {S <: Real}

    @assert class in [:obs, :state, :pseudo] "The allowed class variables are :state, :obs, and :pseudo"

    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    if class == :state
        # States trend should be zero always, so create a DataFrame of zeros
        df           = DataFrame()
        df[!, :date] = quarter_range(date_shockdec_start(m), date_shockdec_end(m))
        for var in keys(m.endogenous_states_augmented)
            df[!, var] .= zero(S)
        end

        return df
    end

    # Otherwise, calculate potentially time-varying trends in DD or DD_pseudo
    hist_regime_inds = regime_indices(m, start_date, end_fcast_date) # do not need to account for ZLB split b/c shocks don't matter for trends
    if hist_regime_inds[1][1] < 1
        hist_regime_inds[1] = 1:hist_regime_inds[1][end]
    end
    if hist_regime_inds[end][end] >= end_index  # if the end index is in the middle of a regime or is past the regime's end
        hist_regime_inds[end] = hist_regime_inds[end][1]:end_index
        fcast_regime_inds = nothing
    else
        fcast_regime_inds = get_fcast_regime_inds(m, forecast_horizons(m; cond_type = cond_type), cond_type,
                                                  start_index = hist_regime_inds[end][end])
        fcast_cutoff = findfirst([regind[end] >= end_index for regind in fcast_regime_inds])
        if isnothing(fcast_cutoff)
            error("The index_shockdec_end(m) occurs past the index of the final forecast date.")
        end
        fcast_regime_inds = fcast_regime_inds[1:fcast_cutoff]
        fcast_regime_inds[end] = fcast_regime_inds[end][1]:end_index
    end
    if length(hist_regime_inds[end]) == 0
        pop!(hist_regime_inds)
    end

    system = compute_system(m;
                            tvis = haskey(get_settings(m), :tvis_information_set))

    trends = Matrix{S}(undef, class == :obs ? n_observables(m) : n_pseudo_observables(m), end_index)
    # Calculate DD or DD_pseudo trend, based on the class
    if class == :obs
        for (reg, inds) in enumerate(hist_regime_inds)
            trends[:, inds] .= system[reg, :DD]
        end

        if !isnothing(fcast_regime_inds)
            for (i, inds) in enumerate(fcast_regime_inds)
                reg = i + length(hist_regime_inds)
                for t in inds
                    trends[:, t] .= system[reg, :DD]
                end
            end
        end
    else
        for (reg, inds) in enumerate(hist_regime_inds)
            trends[:, inds] .= system[reg, :DD_pseudo]
        end

        if !isnothing(fcast_regime_inds)
            for (i, inds) in enumerate(fcast_regime_inds)
                reg = i + length(hist_regime_inds)
                for t in inds
                    trends[:, inds] .= system[reg, :DD_pseudo]
                end
            end
        end
    end

    # Construct DataFrame of trends
    df           = DataFrame()
    df[!, :date] = quarter_range(date_shockdec_start(m), date_shockdec_end(m))
    if start_index != 1
        trends = trends[:, start_index:end]
    end

    if annualize
        for (k, v) in (class == :obs ? m.observables : m.pseudo_observables)
            df[!, k] = 4. .* trends[v, :]
        end
    else
        for (k, v) in (class == :obs ? m.observables : m.pseudo_observables)
            df[!, k] = trends[v, :]
        end
    end

    return df
end

#="""
```
prepare_meansbands_table_irf(mb, var, shock)
```

Returns a `DataFrame` of means and bands for a particular impulse
response function of variable (observable, pseudoobservable, or state)
`v` to shock `s`. Columns are sorted such that the bands are ordered from
smallest to largest, and the means are on the far right. For example,
a MeansBands object containing 50\\% and 68\\% bands would be ordered as
follows: [68\\% lower, 50\\% lower, 50\\% upper, 68\\% upper, mean].

### Inputs
- `mb::MeansBands`: time-series MeansBands object
- `var::Symbol`: an economic variable stored in `mb`. If `mb` stores
  observables, `var` would be an element of `propertynames(m.observables)`. If
  it stores pseudoobservables, `var` would be the name of a
  pseudoobservable defined in the pseudomeasurement equation.
"""
function prepare_meansbands_table_irf(mb::MeansBands, shock::Symbol, var::Symbol)

    @assert get_product(mb) in [:irf] "prepare_meansbands_table_irf can only be used for irfs"
    @assert var in get_vars_means(mb) "$var is not stored in this MeansBands object"

    # get the variable-shock combination we want to print
    # varshock = Symbol["$var" * DSGE_SHOCKDEC_DELIM * "$shock" for var in vars]
    varshock = Symbol("$var" * DSGE_SHOCKDEC_DELIM * "$shock")

    # extract the means and bands for this irf
    df = mb.bands[varshock][map(Symbol, which_density_bands(mb))]
    df[:mean] = mb.means[varshock]

    return df
end

function prepare_meansbands_table_irf(mb::MeansBands, shock::Symbol, vars::Vector{Symbol})

    # Print all vars by default
    if isempty(vars)
        vars = get_variables(mb)
    end

    # Make dictionary to return
    irfs = Dict{Symbol, DataFrame}()
    df = DataFrame()

    # Make tables for each irf
    for var in vars
        @show df
        if isempty(df)
            df = prepare_meansbands_table_irf(mb, shock, var)
        else
            join(df, prepare_meansbands_table_irf(mb, shock, var), on = :date)
        end
    end

    return df
end=#
