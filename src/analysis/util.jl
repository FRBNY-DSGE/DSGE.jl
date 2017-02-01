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
end

function get_product(s::Symbol)
    product = if contains(string(s), "hist")
        :hist
    elseif contains(string(s), "bddforecast4q")
        :bddforecast4q
    elseif contains(string(s), "forecast4q")
        :forecast4q
    elseif contains(string(s), "bddforecast")
        :bddforecast
    elseif contains(string(s), "forecast")
        :forecast
    elseif contains(string(s), "shockdec")
        :shockdec
    elseif contains(string(s), "dettrend")
        :dettrend
    elseif contains(string(s), "trend")
        :trend
    elseif contains(string(s), "irf")
        :irf
    end
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
    if mb.metadata[:product] != :trend && mb.metadata[:product] != :irf
        @printf io "  dates: %s - %s\n" startdate_means(mb) enddate_means(mb)
    end
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

function get_meansbands_input_files{S<:AbstractString}(m::AbstractModel,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     model_string::S = S(""), forecast_string::S = S(""),
                     fileformat = :jld)

    dir = rawpath(m, "forecast", "")
    get_meansbands_input_files(dir, input_type, cond_type, output_vars,
                               model_string = model_string,
                               forecast_string = forecast_string,
                               fileformat = fileformat)
end

function get_meansbands_input_files{S<:AbstractString}(directory::S,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     model_string::S = S(""), forecast_string::S = S(""),
                     fileformat::Symbol = :jld)

    input_files = Dict{Symbol, S}()

    for var in output_vars
        input_files[var] = get_forecast_filename(directory, input_type, cond_type, var,
                                                 model_string = model_string,
                                                 forecast_string = forecast_string,
                                                 fileformat = fileformat)

        if contains(string(var), "4q")
            input_files[var] = replace(input_files[var], "forecast4q", "forecast")
        end
    end

    input_files
end

function get_meansbands_output_files(m::AbstractModel,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     model_string = "", forecast_string = "",
                     fileformat::Symbol = :jld)

    dir = workpath(m, "forecast", "")
    get_meansbands_output_files(dir, input_type, cond_type, output_vars,
                                model_string = model_string,
                                forecast_string = forecast_string, fileformat = fileformat)
end

function get_meansbands_output_files{S<:AbstractString}(directory::S,
                     input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
                     model_string = "", forecast_string = "", fileformat = :jld)

    model_string = S(model_string)
    forecast_string = S(forecast_string)

    mb_output_vars = [symbol("mb$x") for x in output_vars]
    output_files = Dict{Symbol,AbstractString}()

    for var in output_vars
        output_files[var] = get_forecast_filename(directory, input_type, cond_type, symbol("mb$var"),
                                                  model_string = model_string,
                                                  forecast_string = forecast_string,
                                                  fileformat = fileformat)
    end

    output_files
end



"""
```
resize_population_forecast(population_forecast::DataFrame, nperiods::Int;
                                           population_mnemonic::Symbol = Symbol())
```

Extends or shrinks the population forecasts to be `nperiods` in
length. If `population_forecast` must be extended, the last value is
simply repeated as many times as necessary.
"""
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
    for i in 1:length(l1)
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

"""
```
load_population_growth(data_file, forecast_file, population_mnemonic; verbose = :low)
```

Returns `DataFrame`s of growth rates for HP-filtered population data and forecast.
"""
function load_population_growth{S<:AbstractString}(data_file::S, forecast_file::S,
                                                   population_mnemonic::Nullable{Symbol};
                                                   verbose::Symbol = :low)
    if isnull(population_mnemonic)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            warn("No population mnemonic provided")
        end

        return DataFrame(), DataFrame()
    else
        mnemonic = get(population_mnemonic)

        # Read in unfiltered series
        unfiltered_data     = read_population_data(data_file; verbose = :low)
        unfiltered_forecast = read_population_forecast(forecast_file, mnemonic; verbose = :low)

        # HP filter
        data, forecast = transform_population_data(unfiltered_data, unfiltered_forecast,
                                                   mnemonic; verbose = :none)
        dlfiltered_data =
            DataFrame(date = @data(convert(Array{Date}, data[:date])),
                      population_growth = @data(convert(Array{Float64},
                                                        data[:dlfiltered_population_recorded])))
        dlfiltered_forecast =
            DataFrame(date = @data(convert(Array{Date}, forecast[:date])),
                      population_growth = @data(convert(Array{Float64},
                                                        forecast[:dlfiltered_population_forecast])))

        return dlfiltered_data, dlfiltered_forecast
    end
end

"""
```
get_mb_population_series(product, population_mnemonic, population_data, population_forecast, date_list)
```

Returns the appropriate population series for the `product`.
"""
function get_mb_population_series(product::Symbol, population_mnemonic::Nullable{Symbol},
                               population_data::DataFrame, population_forecast::DataFrame,
                               date_list::Vector{Date})
    # Unpack population mnemonic
    mnemonic = isnull(population_mnemonic) ? Symbol() : get(population_mnemonic)

    if product in [:forecast, :bddforecast]

        # For forecasts, we repeat the last forecast period's population
        # forecast until we have n_fcast_periods of population forecasts
        n_fcast_periods = length(date_list)
        population_series = resize_population_forecast(population_forecast, n_fcast_periods,
                                                       population_mnemonic = mnemonic)
        population_series = convert(Vector{Float64}, population_series[mnemonic])

    elseif product in [:shockdec, :dettrend, :trend, :forecast4q, :bddforecast4q]

        if product in [:forecast4q, :bddforecast4q]
            # For forecast4q, we want the last 3 historical periods + the forecast
            # date_list is the date_list for forecast, so date_list[1] corresponds to date_forecast_start.
            start_date = iterate_quarters(date_list[1], -3)
            end_date   = date_list[end]
            start_ind  = find(population_data[:date] .== start_date)[1]
        else
            # For shockdecs, deterministic trend, and trend, we want to
            # make sure population series corresponds with the saved dates.
            start_date = date_list[1]
            end_date   = date_list[end]
            start_ind  = find(population_data[:date] .== start_date)[1]
        end
        population_data = population_data[start_ind:end, mnemonic]

        # Calculate number of periods that are in the future
        n_fcast_periods = if product in [:forecast4q, :bddforecast4q]
            length(date_list)
        else
            length(date_list) - length(population_data)
        end

        # Extend population forecast by the right number of periods
        population_forecast = resize_population_forecast(population_forecast, n_fcast_periods,
                                                         population_mnemonic = mnemonic)
        end_ind = find(population_forecast[:date] .== end_date)[1]

        # Concatenate population histories and forecasts together
        population_series = if isempty(end_ind)
            convert(Vector{Float64}, population_data)
        else
            tmp = [population_data; population_forecast[1:end_ind, mnemonic]]
            convert(Vector{Float64}, tmp)
        end

    elseif product == :hist

        # For history, the population series is just the data
        population_series = convert(Vector{Float64}, population_data[mnemonic])

    end

    return population_series
end

"""
```
get_mb_metadata(input_type, cond_type, output_var, forecast_output_file; forecast_string = "")
```

Returns the `metadata` dictionary from `read_forecast_metadata`, as well as
`mb_metadata`, the dictionary that we will save to the means and bands file.
"""
function get_mb_metadata{S<:AbstractString}(input_type::Symbol, cond_type::Symbol,
                                            output_var::Symbol, forecast_output_file::S;
                                            forecast_string = "")
    class   = get_class(output_var)
    product = get_product(output_var)

    metadata, fcast_output = jldopen(forecast_output_file, "r") do jld
        read_forecast_metadata(jld), read_forecast_output(jld)
    end

    if class == :pseudo
        transforms       = metadata[:pseudoobservable_revtransforms]
        variable_indices = metadata[:pseudoobservable_indices]
    elseif class == :obs
        transforms       = metadata[:observable_revtransforms]
        variable_indices = metadata[:observable_indices]
    elseif class == :shock
        transforms       = metadata[:shock_revtransforms]
        variable_indices = metadata[:shock_indices]
    else
        error("Means and bands are only calculated for observables, pseudo-observables, and shocks")
    end
    date_indices         = product == :irf ? Dict{Date,Int}() : metadata[:date_indices]

    # Make sure date lists are valid. This is vacuously true for trend and IRFs,
    # which are not time-dependent and hence have empty `date_indices`.
    date_list          = collect(keys(date_indices))   # unsorted array of actual dates
    date_indices_order = collect(values(date_indices)) # unsorted array of date indices
    check_consistent_order(date_list, date_indices_order)
    sort!(date_list, by = x -> date_indices[x])

    mb_metadata = Dict{Symbol,Any}(
                   :para            => input_type,
                   :cond_type       => cond_type,
                   :product         => product,
                   :class           => class,
                   :indices         => variable_indices,
                   :forecast_string => forecast_string,
                   :date_inds       => date_indices)

    return fcast_output, metadata, mb_metadata, transforms
end

function get_y0_index(m::AbstractModel, product::Symbol)
    if product in [:forecast, :bddforecast]
        return index_forecast_start(m) - 1
    elseif product in [:forecast4q, :bddforecast4q]
        # We subtract 4 because there is 1 transform that actually
        # needs us to go 4 periods. Later, we can use y0_index + 1
        # to index out the data we need for all the other forecasts.
        return index_forecast_start(m) - 4
    elseif product in [:shockdec, :dettrend, :trend]
        return index_shockdec_start(m) - 1
    elseif product == :hist
        return index_mainsample_start(m) - 1
    elseif product == :irf
        return -1
    else
        error("get_y0_index not implemented for product = $product")
    end
end