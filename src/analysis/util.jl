##########################################
## Extract class and product from a Symbol
##########################################

function get_class(output_var::Symbol)
    s = string(output_var)
    if contains(s, "pseudo")
        :pseudo
    elseif contains(s, "obs")
        :obs
    elseif contains(s, "state")
        :states
    elseif contains(s, "stdshock")
        :stdshocks
    elseif contains(s, "shock")
        :shocks
    else
        error("Invalid output_var: " * s)
    end
end

function get_product(output_var::Symbol)
    s = string(output_var)
    if contains(s, "hist4q")
        :hist4q
    elseif contains(s, "hist")
        :hist
    elseif contains(s, "bddforecast4q")
        :bddforecast4q
    elseif contains(s, "forecast4q")
        :forecast4q
    elseif contains(s, "bddforecastut")
        :bddforecastut
    elseif contains(s, "forecastut")
        :forecastut
    elseif contains(s, "bddforecast")
        :bddforecast
    elseif contains(s, "forecast")
        :forecast
    elseif contains(s, "shockdec")
        :shockdec
    elseif contains(s, "dettrend")
        :dettrend
    elseif contains(s, "trend")
        :trend
    elseif contains(s, "irf")
        :irf
    else
        error("Invalid output_var: " * s)
    end
end

function get_class_longname(class::Symbol)
    if class == :pseudo
        :pseudoobservable
    elseif class == :obs
        :observable
    elseif class == :states
        :state
    elseif class in [:shocks, :stdshocks]
        :shock
    end
end


#####################################
## Preparing data for transformations
#####################################

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
    dates = Base.filter(dr) do x
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
parse_transform(t::Symbol) = eval(Symbol(split(string(t),".")[end]))

"""
```
load_population_growth(data_file, forecast_file, mnemonic;
    use_population_forecast = true, use_hpfilter = true, verbose = :low)
```

Returns `DataFrame`s of growth rates for HP-filtered population data and forecast.
"""
function load_population_growth(data_file::String, forecast_file::String,
                                mnemonic::Symbol;
                                use_population_forecast::Bool = true,
                                use_hpfilter::Bool = true,
                                verbose::Symbol = :low)

    data_verbose = verbose == :none ? :none : :low

    # Read in unfiltered series
    unfiltered_data     = read_population_data(data_file; verbose = data_verbose)
    unfiltered_forecast = if use_population_forecast
        read_population_forecast(forecast_file, mnemonic; verbose = data_verbose)
    else
        DataFrame()
    end

    # HP filter if necessary
    data, forecast = transform_population_data(unfiltered_data, unfiltered_forecast,
                                               mnemonic; use_hpfilter = use_hpfilter,
                                               verbose = :none)
    if use_hpfilter
        data_mnemonic = :dlfiltered_population_recorded
        forecast_mnemonic = :dlfiltered_population_forecast
    else
        data_mnemonic = :dlpopulation_recorded
        forecast_mnemonic = :dlpopulation_forecast
    end

    # Prepare output variables
    data  = data[[:date, data_mnemonic]]
    rename!(data, data_mnemonic, :population_growth)

    if use_population_forecast
        forecast = forecast[[:date, forecast_mnemonic]]
        rename!(forecast, forecast_mnemonic, :population_growth)
    else
        forecast = DataFrame()
    end

    return data, forecast
end

"""
```
get_population_series(mnemonic, population_data, population_forecast,
    start_date, end_date)
```

Returns the population series between `start_date` and `end_date`. If `end_date`
is after the last period in `population_forecast`, the last forecasted
population is repeated until `end_date` using
`resize_population_forecast`. Returns a `Vector{Float64}`.
"""
function get_population_series(mnemonic::Symbol, population_data::DataFrame,
                               population_forecast::DataFrame, start_date::Date,
                               end_date::Date)

    last_historical_date = population_data[end, :date]

    # If no population forecast, use last period of population data instead
    if isempty(population_forecast)
        population_forecast = population_data[end, :]
        last_forecast_date = DSGE.iterate_quarters(last_historical_date, 1)
        population_forecast[1, :date] = last_forecast_date
    else
        last_forecast_date = population_forecast[end, :date]
    end

    # Calculate number of periods that are in the future and extend
    # population forecast by the right number of periods
    if end_date > last_forecast_date
        n_fcast_periods = subtract_quarters(end_date, last_historical_date)
        population_forecast = resize_population_forecast(population_forecast, n_fcast_periods,
                                                         mnemonic)
    end

    population_insample = if start_date <= population_data[end, :date]

        padding = if start_date < population_data[1, :date]
            # Start date is before population data; compute number of NaNs to prepend
            warn("Start date $start_date is before population data begins: prepending NaNs")
            n_nans = subtract_quarters(population_data[1, :date], start_date)

            DataFrame(date = quarter_range(start_date, iterate_quarters(population_data[1,:date], -1)))
        else
            DataFrame()
        end

        unpadded_data = if population_data[1, :date] < end_date < population_data[end, :date]
            # Dates entirely in past
            population_data[start_date .<= population_data[:, :date] .<= end_date, :]
        else
            # Dates span past and forecast
            data  = population_data[start_date .<= population_data[:, :date], :]
            fcast = population_forecast[population_forecast[:, :date] .<= end_date, :]
            vcat(data, fcast)
        end

        padded_data =  vcat(padding, unpadded_data)
        na2nan!(padded_data)
        padded_data

    elseif population_forecast[1, :date] <= start_date <= population_forecast[end, :date]
        # Dates entirely in forecast
        population_forecast[start_date .<= population_forecast[:, :date] .<= end_date, :]
    else
        # start_date comes after population_forecast[end, :date]
        error("Start date $start_date comes after population forecast ends")
    end

    return convert(Vector{Float64}, population_insample[mnemonic])
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

    if product == :irf
        # Return empty vector for IRFs, which don't correspond to real dates
        return Vector{Float64}()
    else
        start_date = if product in [:hist4q, :forecast4q, :bddforecast4q]
            iterate_quarters(date_list[1], -3)
        elseif product in [:hist, :forecast, :bddforecast, :shockdec, :dettrend, :trend]
            date_list[1]
        else
            error("Invalid product: $product")
        end
        end_date = date_list[end]

        return get_population_series(mnemonic, population_data, population_forecast,
                                     start_date, end_date)
    end
end


##############
## Other utils
##############

"""
```
get_mb_metadata(input_type, cond_type, output_var, forecast_output_file; forecast_string = "")
```

Returns the `metadata` dictionary from `read_forecast_metadata`, as well as
`mb_metadata`, the dictionary that we will save to the means and bands file.
"""
function get_mb_metadata(input_type::Symbol, cond_type::Symbol,
                         output_var::Symbol, forecast_output_file::String;
                         forecast_string::String = "")
    class   = get_class(output_var)
    product = get_product(output_var)

    metadata = jldopen(forecast_output_file, "r") do jld
        read_forecast_metadata(jld)
    end

    class_long = get_class_longname(class)
    variable_indices = metadata[Symbol(class_long, "_indices")]
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
    if product in [:forecast, :bddforecast]
        return index_forecast_start(m) - 1
    elseif product in [:forecast4q, :bddforecast4q]
        # We subtract 4 because there is 1 transform that actually
        # needs us to go 4 periods. Later, we can use y0_index + 1
        # to index out the data we need for all the other forecasts.
        return index_forecast_start(m) - 4
    elseif product in [:shockdec, :dettrend, :trend]
        return n_presample_periods(m) + index_shockdec_start(m) - 1
    elseif product in [:hist, :hist4q]
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
