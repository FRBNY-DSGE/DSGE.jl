"""
```
smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, syses::Vector{System},
    kals::Vector{Kalman})

smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System}, kals::Vector{Kalman})
```

Computes and returns the smoothed values of states for every parameter draw.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses::Vector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `kals::Vector{Kalman}`: a vector of `Kalman` objects containing the results of
  the Kalman filter for each draw

### Outputs

- `alpha_hats`: a vector of smoothed states (`alpha_hat`s) returned from the
  smoother specified by `smoother_flag(m)`, one for each system in `syses`
- `eta_hats`: a vector of smoothed shocks (`eta_hat`s) returned from the
  smoother, one for each system in `syses`
"""
function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  df::DataFrame,
                                  syses::Vector{System{S}},
                                  kals::Vector{Kalman{S}})

    data = df_to_matrix(m, df)
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
    
    alpha_hats = [Array(x[1]) for x in out] # to make type stable 
    eta_hats   = [Array(x[2]) for x in out] 

    return alpha_hats, eta_hats
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

