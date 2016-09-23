"""
```
forecast{T<:AbstractFloat}(m::AbstractModel, sys::Vector{System{T}},
    initial_state_draws::Vector{Vector{T}}; shock_distributions::Union{Distribution,
    Matrix{T}} = Matrix{T}())
```

Computes forecasts for all draws, given a model object, system matrices, and a
matrix of shocks or a distribution of shocks

### Inputs

- `m`: model object
- `syses::Vector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `initial_state_draws`: a vector of state vectors in the final historical period
- `shock_distributions`: a `Distribution` to draw shock values from, or
  a matrix specifying the shock innovations in each period

### Outputs

- `states`: 3-dimensional array of size `nstates` x `horizon` x `ndraws`
  consisting of forecasted states for each draw
- `observables`: 3-dimensional array of size `nobs` x `horizon` x `ndraws`
  consisting of forecasted observables for each draw
- `pseudo_observables`: 3-dimensional array of size `npseudo` x `horizon` x `ndraws`
  consisting of forecasted pseudo-observables for each draw
- `shocks`: 3-dimensional array of size `nshocks` x `horizon` x `ndraws`
  consisting of forecasted shocks for each draw
"""
function forecast{T<:AbstractFloat}(m::AbstractModel, syses::Vector{System{T}},
                                    initial_state_draws::Vector{Vector{T}};
                                    shock_distributions::Union{Distribution,
                                    Matrix{T}} = Matrix{T}())

    ndraws = length(syses)

    # Get pseudomeasurement matrices
    Z_pseudo, D_pseudo = if forecast_pseudoobservables(m)
        _, pseudo_mapping = pseudo_measurement(m)
        pseudo_mapping.ZZ, pseudo_mapping.DD
    else
        Matrix{T}(), Vector{T}()
    end
    
    # retrieve settings for forecast
    horizon  = forecast_horizons(m)
    nshocks  = n_shocks_exogenous(m)

    # Unpack everything for call to map/pmap
    TTTs     = [s[:TTT] for s in syses]
    RRRs     = [s[:RRR] for s in syses]
    CCCs     = [s[:CCC] for s in syses]
    ZZs      = [s[:ZZ] for s in syses]
    DDs      = [s[:DD] for s in syses]

    # Prepare copies of these objects due to lack of parallel broadcast functionality
    ZZps     = [Z_pseudo for i in 1:ndraws]
    DDps     = [D_pseudo for i in 1:ndraws]
    horizons = [horizon for i in 1:ndraws]
        
    # set up distribution of shocks if not specified
    # For now, we construct a giant vector of distirbutions of shocks and pass
    # each to compute_forecast.
    #
    # TODO: refactor so that compute_forecast
    # creates its own DegenerateMvNormal based on passing the QQ
    # matrix (which has already been computed/is taking up space)
    # rather than having to copy each Distribution across nodes. This will also be much more
    # space-efficient when forecast_kill_shocks is true.

    shock_distributions = if isempty(shock_distributions)
        if forecast_kill_shocks(m)
            [zeros(nshocks, horizon) for i in 1:ndraws]
        else
            # use t-distributed shocks
            if forecast_tdist_shocks(m)
                [Distributions.TDist(forecast_tdist_df_val(m)) for i in 1:ndraws]
            # use normally distributed shocks
            else
                shock_distributions = Vector{DSGE.DegenerateMvNormal}(ndraws)
                for i = 1:ndraws
                    shock_distributions[i] = DSGE.DegenerateMvNormal(zeros(nshocks),sqrt(syses[i][:QQ]))
                end
                shock_distributions
            end
        end
    end

    if use_parallel_workers(m)
        mapfcn = pmap
    else
        mapfcn = map
    end

    # Go to work!
    forecasts =  mapfcn(DSGE.compute_forecast, TTTs, RRRs, CCCs, ZZs, DDs,
                        horizons, shock_distributions, initial_state_draws,
                        ZZps, DDps)

    # Unpack returned vector of tuples
    states             = [forecast[1]::Matrix{T} for forecast in forecasts]
    observables        = [forecast[2]::Matrix{T} for forecast in forecasts]
    pseudo_observables = [forecast[3]::Matrix{T} for forecast in forecasts]
    shocks             = [forecast[4]::Matrix{T} for forecast in forecasts]

    # Splat vectors of matrices into 3-D arrays
    states             = cat(3, states...)
    observables        = cat(3, observables...)
    pseudo_observables = cat(3, pseudo_observables...)
    shocks             = cat(3, shocks...)

    return states, observables, pseudo_observables, shocks
end

"""
```
compute_forecast(T, R, C, Z, D, forecast_horizons,
    shocks, z,  Z_pseudo, D_pseudo)
```

### Inputs

- `T`, `R`, `C`: transition equation matrices
- `Z`, `D`: observation equation matrices
- `Z_pseudo`, `D_pseudo`: matrices mapping states to pseudo-observables
- `forecast_horizons`: number of quarters ahead to forecast output
- `shocks`: joint distribution (type `Distribution`) from which to draw
  time-invariant shocks or matrix of drawn shocks (size `nshocks` x
  `forecast_horizons`)
- `z`: state vector at time `T`, i.e. at the beginning of the forecast

### Outputs

`compute_forecast` returns a 4-tuple of forecast outputs, whose elements have
sizes:

- `states`: `nstates` x `forecast_horizons`
- `observables`: `nobs` x `forecast_horizons`
- `pseudo_observables`: `npseudo` x `forecast_horizons`
- `shocks`: `nshocks` x `forecast_horizons`
"""
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},
                                            forecast_horizons::Int,
                                            shocks::Matrix{S},
                                            z::Vector{S},
                                            Z_pseudo::Matrix{S}=Matrix{S}(),
                                            D_pseudo::Vector{S}=Vector{S}())

    if forecast_horizons <= 0
        throw(DomainError())
    end
                                    
    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobservables = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    states       = zeros(nstates, forecast_horizons)
    
    # Define our iteration function
    iterate(z_t1, ϵ_t) = C + T*z_t1 + R*ϵ_t

    # Iterate first period
    states[:, 1] = iterate(z, shocks[:, 1])
    
    # Iterate remaining periods
    for t in 2:forecast_horizons
        states[:, t] = iterate(states[:, t-1], shocks[:, t])
    end

    # Apply observation and pseudo-observation equations
    observables        = D        .+ Z        * states
    pseudo_observables = if isempty(Z_pseudo) || isempty(D_pseudo)
        Matrix{S}()
    else
        D_pseudo .+ Z_pseudo * states
    end
    
    # Return a dictionary of forecasts
    return states, observables, pseudo_observables, shocks
end

# Utility method to actually draw shocks
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},  
                                            forecast_horizons::Int,
                                            dist::Distribution,
                                            z::Vector{S},
					    Z_pseudo::Matrix{S}=Matrix{S}(),
                                            D_pseudo::Vector{S}=Vector{S}())

    if forecast_horizons <= 0
        throw(DomainError())
    end

    nshocks = size(R, 2)
    shocks = zeros(nshocks, forecast_horizons)

    for t in 1:forecast_horizons
        shocks[:, t] = rand(dist)
    end

    compute_forecast(T, R, C, Z, D, forecast_horizons, shocks, z,
                     Z_pseudo, D_pseudo)
end

# I'm imagining that a Forecast object could be returned from
# compute_forecast rather than a dictionary. It could look something like the
# type outlined below for a single draw.
# Perhaps we could also add some fancy indexing to be able to index by names of states/observables/etc.
# Question becomes: where do we store the list of observables? In the
# model object? In a separate forecastSettings vector that we also
# pass to forecast?
# immutable Forecast{T<:AbstractFloat}
#     states::Matrix{T}
#     observables::Matrix{T}
#     pseudoobservables::Matrix{T}
#     shocks::Matrix{T}
# end

