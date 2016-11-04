"""
```
compute_means_bands_shockdec(fcast_output, transforms, var_inds, shock_inds, date_list,
                             data = [], population_forecast = DataFrame(), hist_end_index = 0)
```


### Inputs

- `fcast_output`: an `ndraws` x `nvars` x `nperiods` x `nshocks`
  array of shock decompositions from the output of the forecast

"""
function compute_means_bands_shockdec{T<:AbstractFloat}(fcast_output::Array{T},
                                                        transforms::Dict{Symbol,Symbol},
                                                        variable_indices::Dict{Symbol,Int},
                                                        shock_inds::Dict{Symbol,Int},
                                                        date_list::Vector{Date};
                                                        data::Matrix{T} = Matrix{T}(),
                                                        population_series = Array{T},
                                                        hist_end_index::Int = 0,
                                                        density_bands::Array{Float64} = [0.5,0.6,0.7,0.8,0.9])


    # set up means and bands structures
    sort!(date_list)
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()

    # for each element of shock x variable (variable = each pseudoobs, obs, or state):
    # 1. apply the appropriate transform
    # 2. add to DataFrame
    for (shock, shock_ind) in shock_inds
        for (var, var_ind) in variable_indices
            transform = DSGE.parse_transform(transforms[var])
            ex = if transform in [:logtopct_annualized_percapita]
                Expr(:call, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2), population_series)
            elseif transform in [:loglevelto4qpct_annualized_percapita]
                Expr(:call, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2),
                     data[ind,:], hist_end_index, population_series)
            else
                Expr(:call, :map, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2))
            end

            transformed_fcast_output = eval(ex)

            means[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = vec(mean(transformed_fcast_output,1))
            bands_one = find_density_bands(transformed_fcast_output, density_bands, minimize=false)
            bands_one[:date] = date_list
            bands[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = bands_one
        end
    end

    return means, bands
end


function get_shockdec_means(mb::MeansBands,shock::Symbol; vars::Vector{Symbol}=Vector{Symbol}())

    vars = if isempty(vars)
        collect(names(mb.means))[find([contains(string(col), string(shock)) for col in names(mb.means)])]
    else
        vars
    end

    out = DataFrame()
    for var in vars
        varname = split(string(var), DSGE_SHOCKDEC_DELIM)[1]
        out[symbol(varname)] = mb.means[var]
    end

    out
end


function get_shockdec_bands(mb::MeansBands, shock::Symbol; vars::Vector{Symbol}=Vector{Symbol}(), bands=[])

    # If var is not supplied, return all variables
    if isempty(vars)
        vars = collect(keys(mb.bands))[find([contains(string(col), string(shock)) for col in keys(mb.bands)])]
    end

    bands_keys = if isempty(bands)
        names(mb.bands[vars[1]])
    else
        [[symbol("$(100x)% LB") for x in bands]; [symbol("$(100x)% UB") for x in bands]]
    end

    out = Dict{Symbol, DataFrame}()
    for var in vars
        out[var] = mb.bands[var][bands_keys]
    end

    out
end
