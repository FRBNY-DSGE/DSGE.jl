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
                                  forecast_settings::Vector{Setting},
                                  data::Matrix{AbstractFloat},
                                  sys::Vector{System},
                                  kal::Vector{Kalman})

    # numbers of useful things
    ndraws   = n_draws(m)
    @assert length(sys) == ndraws
    @assert length(kal) == ndraws

    # Extract results from Kalman filter: need pred, vpred

    # Extract settings from model: peachcount, psize, n_anticipated_shocks, n_anticipated_lags
    
    # Make output vector
    smoothed = Vector{KalmanSmooth}(ndraws)
    
    # Parallelize
    # Think about whether level 1 functions should return KalmanSmooth
    # objects and level 2 should return a vector of them.
    for i = 1:ndraws
        if disturbance_smoother_flag
            smoothed[i] = disturbance_smoother(...)
        else
            smoothed[i] = kalman_smoother(...)
        end
    end

    return smoothed
end


