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

- `smoothed_states`: a vector of `alpha_hat`s returned from the smoother
  specified by `smoother_flag(m)`, one for each system in `syses`
- `smoothed_shocks`: a vector of `eta_hat`s returned from the smoother, one for
  each system in `syses`
"""
function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  df::DataFrame,
                                  syses::Vector{System},
                                  kals::Vector{Kalman})
    data = df_to_matrix(m, df)'
    smooth(m, data, syses, kals)
end

function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S}
                                  syses::Vector{System},
                                  kals::Vector{Kalman})

    # numbers of useful things
    ndraws = length(syses)
    @assert length(kals) == ndraws

    # Broadcast models and data matrices 
    models = fill(m, ndraws)
    datas = fill(data, ndraws)
    include_presamples = if include_presample
        fill(IncludePresample(), ndraws)
    else
        fill(ExcludePresample(), ndraws)
    end

    # Call smooth over all draws
    if use_parallel_workers(m) && nworkers() > 1
        println("Using pmap")
        mapfcn = pmap
    else
        mapfcn = map
    end    
    out = mapfcn(DSGE.tricky_smooth, include_presamples, models, datas, syses, kals)

    smoothed_states = [Array(x[1]) for x in out]  # to make type stable
    smoothed_shocks = [Array(x[2]) for x in out]  
    
    return smoothed_states, smoothed_shocks
end

tricky_smooth(::IncludePresample, m::AbstractModel, data::Matrix, sys::System, kal::Kalman) = 
    smooth(m, data, sys, kal; include_presample = true)
tricky_smooth(::ExcludePresample, m::AbstractModel, data::Matrix, sys::System, kal::Kalman) = 
    smooth(m, data, sys, kal; include_presample = false)

function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{AbstractFloat},
                                  sys::System,
                                  kal::Kalman)

    alpha_hat, eta_hat = if smoother_flag(m) == :kalman
        kalman_smoother(m, data, sys, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif smoother_flag(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, sys, kal[:z0], kal[:vz0])
    end

    return alpha_hat, eta_hat
end

