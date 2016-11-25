"""
```
compute_means_bands_shockdec(fcast_output, transforms, var_inds, shock_inds, date_list,
                             data = [], population_forecast = DataFrame(), y0_index = 0)
```


### Inputs

* `fcast_output`: an `ndraws` x `nvars` x `nperiods` x `nshocks`
  array of shock decompositions from the output of the forecast

* `data`: an `nperiods` x `nobs` matrix of data
"""
function compute_means_bands_shockdec{T<:AbstractFloat}(fcast_output::Array{T},
                                                        transforms::Dict{Symbol,Symbol},
                                                        variable_indices::Dict{Symbol,Int},
                                                        shock_inds::Dict{Symbol,Int},
                                                        date_list::Vector{Date};
                                                        data::Matrix{T} = Matrix{T}(),
                                                        population_series = Vector{T}(),
                                                        y0_index::Nullable{Int} = Nullable{Int}(),
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
                     data[var_ind, get(y0_index)], population_series)
            else
                Expr(:call, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2))
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
