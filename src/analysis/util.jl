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
    elseif contains(string(s), "bddforecastq4q4")
        :bddforecastq4q4
    elseif contains(string(s), "forecastq4q4")
        :forecastq4q4
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

function get_class_longname(class::Symbol)
    longname = if class == :pseudo
        :pseudoobservable
    elseif class == :obs
        :observable
    elseif class == :state
        :state
    elseif class == :shock
        :shock
    end
end


#########################################
## Preparing data for transformations
#########################################
"""
```
resize_population_forecast(population_forecast::DataFrame, nperiods::Int,
    mnemonic::Symbol)
```

Extends or shrinks the population forecasts to be `nperiods` in
length. If `population_forecast` must be extended, the last value is
simply repeated as many times as necessary.
"""
function resize_population_forecast(population_forecast::DataFrame, nperiods::Int,
                                    mnemonic::Symbol)

    # number of periods to extend population forecast
    n_filler_periods = nperiods - size(population_forecast,1)

    # create date range. There are on average 91.25 days in a quarter.
    last_provided = population_forecast[end,:date]
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
load_population_growth(data_file, forecast_file, mnemonic; verbose = :low)
```

Returns `DataFrame`s of growth rates for HP-filtered population data and forecast.
"""
function load_population_growth{S<:AbstractString}(data_file::S, forecast_file::S,
                                                   mnemonic::Symbol;
                                                   verbose::Symbol = :low)
    data_verbose = verbose == :none ? :none : :low

    # Read in unfiltered series
    unfiltered_data     = read_population_data(data_file; verbose = data_verbose)
    unfiltered_forecast = read_population_forecast(forecast_file, mnemonic; verbose = data_verbose)

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

"""
```
get_mb_population_series(product, mnemonic, population_data, population_forecast, date_list)
```

Returns the appropriate population series for the `product`.
"""
function get_mb_population_series(product::Symbol, mnemonic::Symbol,
                                  population_data::DataFrame, population_forecast::DataFrame,
                                  date_list::Vector{Date})

    if product in [:forecast, :bddforecast]

        # For forecasts, we repeat the last forecast period's population
        # forecast until we have n_fcast_periods of population forecasts
        n_fcast_periods = length(date_list)
        population_series = resize_population_forecast(population_forecast, n_fcast_periods,
                                                       mnemonic)
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
                                                         mnemonic)
        end_ind = if isempty(population_forecast[:date])
            []
        else
            find(population_forecast[:date] .== end_date)[1]
        end

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

    elseif product == :irf

        # Return empty vector for IRFs, which don't correspond to real dates
        population_series = Vector{Float64}()

    else

        error("Invalid product: $product")

    end

    return population_series
end


#########################################
## Other utils
#########################################

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

    metadata = jldopen(forecast_output_file, "r") do jld
        read_forecast_metadata(jld)
    end

    class_long = get_class_longname(class)
    variable_indices = metadata[symbol("$(class_long)_indices")]
    date_indices     = product == :irf ? Dict{Date,Int}() : metadata[:date_indices]

    # Make sure date lists are valid. This is vacuously true for and IRFs, which
    # are not time-dependent and hence have empty `date_indices`.
    date_list          = collect(keys(date_indices))   # unsorted array of actual dates
    date_indices_order = collect(values(date_indices)) # unsorted array of date indices
    check_consistent_order(date_list, date_indices_order)

    mb_metadata = Dict{Symbol,Any}(
                   :para            => input_type,
                   :cond_type       => cond_type,
                   :product         => product,
                   :class           => class,
                   :indices         => sort(variable_indices, by = x -> variable_indices[x]),
                   :date_inds       => sort(date_indices, by = x -> date_indices[x]),
                   :forecast_string => forecast_string)

    return metadata, mb_metadata
end

"""
```
get_y0_index(m::AbstractModel, product::Symbol)
```

Returns the index of the period *before* the start of saved results
for a particular forecast-step product. For regular forecasts, this is
the last period for which we have full-quarter data. For running 4q
forecasts, this corresponds to 4 periods before the first forecast
period. For shockdecs, this is one period before
`index_shockdec_start(m)`, and for histories, it is the last period of
the presample.
"""
function get_y0_index(m::AbstractModel, product::Symbol)
    if product in [:forecast, :bddforecast, :forecast4q, :bddforecast4q, :forecastq4q4, :bddforecastq4q4]
        return index_forecast_start(m) - 1
    elseif product in [:forecast4q, :bddforecast4q]
        # We subtract 4 because there is 1 transform that actually
        # needs us to go 4 periods. Later, we can use y0_index + 1
        # to index out the data we need for all the other forecasts.
        return index_forecast_start(m) - 4
    elseif product in [:shockdec, :dettrend, :trend]
        return n_presample_periods(m) + index_shockdec_start(m) - 1
    elseif product == :hist
        return index_mainsample_start(m) - 1
    elseif product == :irf
        return -1
    else
        error("get_y0_index not implemented for product = $product")
    end
end


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

