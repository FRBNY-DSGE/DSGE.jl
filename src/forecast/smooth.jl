"""
```
smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, syses::Vector{System},
    kals::Vector{Kalman}; cond_type::Symbol = :none)

smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System}, kals::Vector{Kalman})
```

Computes and returns the smoothed values of states and shocks for every
parameter draw.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses::Vector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `kals::Vector{Kalman}`: a vector of `Kalman` objects containing the results of
  the Kalman filter for each draw

### Outputs

- `states`: 3-dimensional array of size `nstates` x `hist_periods` x `ndraws`
  consisting of smoothed states for each draw
- `shocks`: 3-dimensional array of size `nshocks` x `hist_periods` x `ndraws`
  consisting of smoothed shocks for each draw

where `states` and `shocks` are returned from the smoother specified by
`smoother_flag(m)`.
"""
function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  df::DataFrame,
                                  syses::Vector{System{S}},
                                  kals::Vector{Kalman{S}};
                                  cond_type::Symbol = :none)

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth(m, data, syses, kals)
end

function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  syses::Vector{System{S}},
                                  kals::Vector{Kalman{S}})

    # numbers of useful things
    ndraws = length(syses)
    @assert length(kals) == ndraws

    # Broadcast models and data matrices 
    models = fill(m, ndraws)
    datas = fill(data, ndraws)

    # Call smooth over all draws
    if use_parallel_workers(m) && nworkers() > 1
        mapfcn = pmap
    else
        mapfcn = map
    end    
    out = mapfcn(DSGE.smooth, models, datas, syses, kals)

    # Unpack returned vector of tuples
    states = [x[1]::Matrix{S} for x in out]
    shocks = [x[2]::Matrix{S} for x in out]

    # Splat vectors of matrices into 3-D arrays
    states = cat(3, states...)
    shocks = cat(3, shocks...)

    return states, shocks
end

function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  sys::System{S},
                                  kal::Kalman{S})

    alpha_hat, eta_hat = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, sys, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, sys, kal[:z0], kal[:vz0])
    end

    return alpha_hat, eta_hat
end

