
"""
```
compute_means_bands_irf(fcast_output, shock_inds)
```


### Inputs

* `fcast_output`: an `ndraws` x `nvars` x `nperiods` x `nshocks`
  array of shock decompositions from the output of the forecast
* `shock_inds`: a dictionary mapping shocks to indices
* `variables_inds`: a dictionary mapping variables to indices
"""
function compute_means_bands_irf{T<:AbstractFloat}(fcast_output::Array{T}, 
                                                  shock_inds::Dict{Symbol,Int},
                                                  variable_inds::Dict{Symbol,Int};
                                                  density_bands::Array{Float64} = [0.5,0.6,0.7,0.8,0.9])
    # draws x states x periods x shocks
    means = DataFrame()
    bands = Dict{Symbol,DataFrame}()

    # for each element of shock x variable, calculate the mean and bands
    for (shock, shock_ind) in shock_inds
        for (var, var_ind) in variable_inds
            fcast_series = squeeze(fcast_output[:,var_ind,:,shock_ind], 2)
            means[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = vec(mean(fcast_series,1))
            bands_one = find_density_bands(fcast_series, density_bands, minimize=false)
            bands[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = bands_one
        end
    end

    return means, bands
end