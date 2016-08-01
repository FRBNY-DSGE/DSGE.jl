"""
```
smooth(...)
```
Computes and returns the smoothed values of states for every parameter draw.

Inputs
------

- `m`: model object
- `data`: matrix of data for observables
- `sys::Vector{System}`: a vector of `System` objects specifying state-space system matrices for each draw
- `kal::Vector{Kalman}`: a vector of `Kalman` objects containing the results of the Kalman filter for each draw

Outputs
-------
- A vector of `KalmanSmooth` objects containing the smoothed states and shocks for each draw
"""
function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{AbstractFloat},
                                  syses::Vector{System},
                                  kals::Vector{Kalman})

    # numbers of useful things
    ndraws = length(syses)
    @assert length(kals) == ndraws

    # Broadcast models and data matrices 
    models = fill(m, ndraws)
    datas = fill(data, ndraws)

    # Call smooth over all draws
    if use_parallel_workers(m) && nworkers() > 1
        println("Using pmap")
        mapfcn = pmap
    else
        mapfcn = map
    end    
    out = mapfcn(smooth, models, datas, syses)

    smoothed_states = [Array(x[1]) for x in out]  # to make type stable
    smoothed_shocks = [Array(x[2]) for x in out]  
    
    return smoothed_states, smoothed_shocks
end

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

