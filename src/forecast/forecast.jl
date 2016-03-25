"""
`forecast(...)`

Computes forecasts for all draws, given data, system matrices, matrix
of shocks or distribution of shocks, model object

Inputs
------
- `m`: model object
- `data`: matrix of data for observables
- `sys::Vector{System}`: a vector of `System` objects specifying
  state-space system matrices for each draw
- `z_T`: state vector in final historical period
- `shockdistribution`: a `Distribution` to draw shock values from, or
  a matrix specifying the shock innovations in each period
- `observables::Vector{Symbol}`: a list of observables to forecast

Outputs
-------
- A vector of `Forecast` objects, which contain the forecasted values of
states, observables, and pseudoobserables, as well as the values of
shock innovations.

"""
function forecast(m::AbstractModel,
                  data::Matrix{AbstractFloat},
                  sys::Vector{System},    
                  z_T::Vector{AbstractFloat}, 
                  shockdistribution::Distribution,  
                  observables::Vector{Symbol})
    
    # numbers of useful things
    ndraws   = n_draws(m)
    @assert length(sys) == ndraws
    
    # make matrix of forecasts and forecast shock values
    forecasts = Vector{Forecast}(ndraws)
    
    # Parallelize over all draws
    for i = 1:ndraws
        forecasts[i] = compute_forecast(...)
    end

    return forecasts
end

# I'm imagining that a Forecast object looks like the following for a single draw.
# Perhaps we could also add some fancy indexing to be able to index by names of states/observables/etc.
# Question becomes: where do we store the list of observables? In the
# model object? In a separate forecastSettings vector that we also
# pass to forecast?
immutable Forecast{T<:AbstractFloat}
    states::Matrix{T}
    observables::Matrix{T}
    pseudoobservables::Matrix{T}
    shocks::Matrix{T}
end
