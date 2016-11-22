#########################################
## Extract class and product from a symbol
#########################################

function get_class(s::Symbol)
    class = if contains(string(s), "pseudo")
        :pseudo
    elseif contains(string(s), "obs")
        :obs
    elseif contains(string(s), "state")
        :state
    elseif contains(string(s), "shock")
        :shock
    end

    class
end

function get_product(s::Symbol)
    product = if contains(string(s), "hist")
        :hist
    elseif contains(string(s), "forecast")
        :forecast
    elseif contains(string(s), "shockdec")
        :shockdec
    elseif contains(string(s), "dettrend")
        :dettrend
    elseif contains(string(s), "trend")
        :trend
    end

    product
end



#########################################
## Useful methods for MeansBands objects
#########################################


class(mb::MeansBands) = mb.metadata[:class]
product(mb::MeansBands) = mb.metadata[:product]
cond_type(mb::MeansBands) = mb.metadata[:cond_type]
para(mb::MeansBands) = mb.metadata[:para]

function Base.show(io::IO, mb::MeansBands)
    @printf io "MeansBands\n"
    @printf io "  class: %s\n" class(mb)
    @printf io "  product: %s\n" product(mb)
    @printf io "  cond: %s\n" cond_type(mb)
    @printf io "  para: %s\n" para(mb)
    @printf io "  dates: %s - %s\n" startdate_means(mb) enddate_means(mb)
    @printf io "  # of variables: %s\n" n_vars_means(mb)
    @printf io "  bands: %s\n" get_density_bands(mb, uniqueify=true)
end

"""
```
function read_mb{S<:AbstractString}(fn::S)
```
Read in a `MeansBands` object saved in `fn`.
"""
function read_mb{S<:AbstractString}(fn::S)
    @assert isfile(fn) "File $fn could not be found"
    mb = jldopen(fn, "r") do f
        read(f, "mb")
    end
    mb
end

###################################
## MEANS
###################################

"""
```
n_vars_means(mb::MeansBands)
````

Get number of variables (`:y_t`, `:OutputGap`, etc) in `mb.means`
"""
function n_vars_means(mb::MeansBands)
    length(get_vars_means(mb))
end

"""
```
get_vars_means(mb::MeansBands)
````

Get variables (`:y_t`, `:OutputGap`, etc) in `mb.means`, sorted by `mb.metadata[:indices]`
"""
function get_vars_means(mb::MeansBands)
    unzip(sort_by_value(mb.metadata[:indices]))[2]
end


"""
```
n_periods_means(mb::MeansBands)
```

Get number of periods in `mb.means`
"""
n_periods_means(mb::MeansBands) = size(mb.means,1)

"""
```
startdate_means(mb::MeansBands)
```

Get first period in`mb.means`. Assumes `mb.means[product]` is already sorted by date.
"""
startdate_means(mb::MeansBands) = mb.means[:date][1]

"""
```
enddate_means(mb::MeansBands)
```

Get last period for which `mb` stores means. Assumes `mb.means[product]` is already sorted by date.
"""
enddate_means(mb::MeansBands) = mb.means[:date][end]


"""
```
get_shockdec_means(mb::MeansBands, var::Symbol; shocks::Vector{Symbol}=Vector{Symbol}())
```

Return the mean value of each shock requested in the shock decomposition of a particular variable.
If `shocks` is empty, returns all shocks.
"""

function get_shockdec_means(mb::MeansBands, var::Symbol; shocks::Vector{Symbol}=Vector{Symbol}())

    # Extract the subset of columns relating to the variable `var` and the shocks listed in `shocks.`
    # If `shocks` not provided, give all the shocks
    var_cols = collect(names(mb.means))[find([contains(string(col), string(var)) for col in names(mb.means)])]
    if !isempty(shocks)
        var_cols = [col -> contains(string(col), string(shock)) ? col : nothing for shock in shocks]
    end

    # Make a new DataFrame with column the column names
    out = DataFrame()
    for col in var_cols
        shockname = split(string(col), DSGE_SHOCKDEC_DELIM)[2]
        out[symbol(shockname)] = mb.means[col]
    end

    out
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
`product` (`hist`, `forecast`, `shockdec`, etc).
"""
function n_periods_bands(mb::MeansBands)
    size(mb.bands[collect(keys(mb.bands))[1]],1)
end

"""
```
startdate_bands(mb::MeansBands)
```

Get first period for which `mb` stores bands. Assumes `mb.bands` is already sorted by date.
"""
startdate_bands(mb::MeansBands) = mb.bands[collect(keys(mb.bands))][:date][1]

"""
```
enddate_bands(mb::MeansBands)
```

Get last period in `mb.bands`. Assumes `mb.bands` is already sorted by date.
"""
enddate_bands(mb::MeansBands) = mb.bands[:date]

"""
```
get_density_bands(mb, uniqueify=false)
```

Return a list of the bands stored in mb.bands. If `uniqueify=true`,
strips "upper" and "lower" band tags and returns unique list of percentage values.
"""
function get_density_bands(mb::MeansBands; uniqueify=false, ordered=true)

    # extract one of the keys in mb.bands
    var  = collect(keys(mb.bands))[1]

    # get all the columns in the corresponding dataframe that aren't dates
    strs = map(string,names(mb.bands[var]))
    strs = setdiff(strs, ["date"])

    lowers = strs[map(ismatch, repmat([r"LB"], length(strs)), strs)]
    uppers = strs[map(ismatch, repmat([r"UB"], length(strs)), strs)]

    # sort
    if ordered
        sort!(lowers, rev=true)
        sort!(uppers)
    end

    # return both upper and lower bands, or just percents, as desired
    strs = if uniqueify
        sort(unique([split(x, " ")[1] for x in [lowers; uppers]]))
    else
        [lowers; uppers]
    end

    return strs
end


"""
```
get_shockdec_bands(mb::MeansBands, var::Symbol;
       shocks::Vector{Symbol}=Vector{Symbol}(), bands::Vector{Symbol}()=Vector{Symbol}())
```

Return bands stored in this MeansBands object for all variables.

### Inputs
- `mb`: MeansBands object
- `var`: the variable of interest (eg the state `:y_t`, or observable `:obs_hours`)

### Optional arguments
- `shocks`: subset of shock names for which to return bands. If empty, `get_shockdec_bands` returns all bands.
- `bands`: subset of bands stored in the DataFrames of `mb.bands` to return.

### Outputs

A `Dict{Symbol, DataFrame}` mapping names of shocks to the bands of `var` corresponding to each shock.
"""
function get_shockdec_bands(mb::MeansBands, var::Symbol;
                            shocks::Vector{Symbol}=Vector{Symbol}(), bands::Vector{Symbol}=Vector{Symbol}())

    # Extract the subset of columns relating to the variable `var` and the shocks listed in `shocks.`
    # If `shocks` not provided, give all the shocks
    var_cols = collect(keys(mb.bands))[find([contains(string(col), string(var)) for col in keys(mb.bands)])]
    if !isempty(shocks)
        var_cols = [col -> contains(string(col), string(shock)) ? col : nothing for shock in shocks]
    end

    # Extract the subset of bands we want to return. Return all bands if `bands` not provided.
    bands_keys = if isempty(bands)
        names(mb.bands[var_cols[1]])
    else
        [[symbol("$(100x)% LB") for x in bands]; [symbol("$(100x)% UB") for x in bands]]
    end

    # Make a new dictionary mapping shock names to bands
    out = Dict{Symbol, DataFrame}()
    for col in var_cols
        shockname = split(string(col), DSGE_SHOCKDEC_DELIM)[2]
        out[shockname] = mb.bands[col][bands_keys]
    end

    out
end


#####################################
## OTHER UTILS
#####################################

function resize_population_forecast(population_forecast::DataFrame, nperiods::Int;
                                           population_mnemonic::Symbol = Symbol())

    # number of periods to extend population forecast
    n_filler_periods = nperiods - size(population_forecast,1)

    # Extract population mnemonic from Nullable object (if not null). If null,
    # take a guess or throw an error if you can't tell.
    mnemonic = if population_mnemonic == Symbol()
        if size(population_forecast,2) > 2
            error("Please indicate which column contains population
                forecasts using the population_mnemonic keyword argument")
        end

        setdiff(names(population_forecast), [:date])[1]
    else
        population_mnemonic
    end

    last_provided = population_forecast[end,:date]

    # create date range. There are on average 91.25 days in a quarter.
    dr = last_provided:(last_provided+Dates.Day(93 * n_filler_periods))

    islastdayofquarter = x->Dates.lastdayofquarter(x) == x
    dates = recur(dr) do x
        islastdayofquarter(x)
    end

    # first element of dates is the last quarter of population forecasts supplied, so we
    # cut it off here
    dates = dates[2:n_filler_periods+1]

    extra = DataFrame()
    extra[:date] = dates

    # resize population forecast by adding n_filler_periods to the given forecast.
    resized = if n_filler_periods > 0

        extra_stuff = fill(population_forecast[end,mnemonic], n_filler_periods)
        extra[mnemonic] = extra_stuff

        [population_forecast; extra]
    elseif n_filler_periods < 0
        population_forecast[1:nperiods,:]
    else
        population_forecast
    end

    resized
end

"""
```
parse_transform(t::Symbol)
```

Parse the module name out of a Symbol to recover the transform associated with
an observable or pseudoobservable. Returns a function.
"""
parse_transform(t::Symbol) = eval(symbol(split(string(t),".")[end]))

"""
```
check_consistent_order(l1, l2)
```

Checks to make sure that l1 and l2 are ordered consistently.
"""
function check_consistent_order(l1, l2)
    @assert length(l1) == length(l2)

    # cache old pairs in a Dict
    original_pairs = Dict{eltype(l1), eltype(l2)}()
    for (i,item) in enumerate(l1)
        original_pairs[item] = l2[i]
    end

    # sort
    l1_sorted = sort(l1)
    l2_sorted = sort(l2)

    # make sure each pair has same pair as before
    for i in length(l1)
        @assert original_pairs[l1_sorted[i]] == l2_sorted[i] "Lists not consistently ordered"
    end

    return true
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
    map(symbol, split(string(s), DSGE_SHOCKDEC_DELIM))
end

"""
```
sort_by_value(d::Dict)
```

Sort a dictionary by value.
"""
function sort_by_value(d::Dict)
    sort(collect(zip(values(d),keys(d))))
end

function unzip{T<:Tuple}(A::Array{T})
    res = map(x -> x[], T.parameters)
    res_len = length(res)
    for t in A
        for i in 1:res_len
            push!(res[i], t[i])
        end
    end
    res
end

