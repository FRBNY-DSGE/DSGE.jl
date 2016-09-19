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

-`states`: vector of length `ndraws`, whose elements are the `nstates` x
 `horizon` matrices of forecasted states
-`observables`: vector of length `ndraws`, whose elements are the `nobs` x
 `horizon` matrices of forecasted observables
-`pseudo_observables`: vector of length `ndraws`, whose elements are the
 `npseudo` x `horizon` matrices of forecasted pseudo-observables
-`shocks`: vector of length `ndraws`, whose elements are the `npseudo` x
 `horizon` matrices of shock innovations
"""
function forecast{T<:AbstractFloat}(m::AbstractModel, syses::Vector{System{T}},
                                    initial_state_draws::Vector{Vector{T}};
                                    shock_distributions::Union{Distribution,
                                    Matrix{T}} = Matrix{T}())

    ndraws = length(syses)

    # for now, we are ignoring pseudo-observables so these can be empty
    Z_pseudo = zeros(T, 12, n_states_augmented(m))
    D_pseudo = zeros(T, 12)

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

    # Forecast the states
    if use_parallel_workers(m)
        mapfcn = pmap
    else
        mapfcn = map
    end

    # Go to work!
    forecasts =  mapfcn(DSGE.compute_forecast, TTTs, RRRs, CCCs, ZZs, DDs, ZZps, DDps,
                              horizons, shock_distributions, initial_state_draws)

    # unpack the giant vector of dictionaries that gets returned
    states      = [x[:states]::Matrix{T} for x in forecasts]
    observables = [x[:observables]::Matrix{T} for x in forecasts]
    pseudo      = [x[:pseudo_observables]::Matrix{T} for x in forecasts]
    shocks      = [x[:shocks]::Matrix{T} for x in forecasts]

    return states, observables, pseudo, shocks
end

"""
```
compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
    shocks, z)
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

`compute_forecast` returns a dictionary of forecast outputs, with keys:

-`:states`
-`:observables`
-`:pseudo_observables`
-`:shocks`
"""
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},
                                            Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
                                            forecast_horizons::Int,
                                            shocks::Matrix{S},
                                            z::Vector{S})

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
    pseudo_observables = D_pseudo .+ Z_pseudo * states
    
    # Return a dictionary of forecasts
    Dict{Symbol, Matrix{S}}(
        :states             => states,
        :observables        => observables,
        :pseudo_observables => pseudo_observables,
        :shocks             => shocks)
end

# Utility method to actually draw shocks
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},          
                                            Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
                                            forecast_horizons::Int,
                                            dist::Distribution,
                                            z::Vector{S})

    if forecast_horizons <= 0
        throw(DomainError())
    end

    nshocks = size(R, 2)
    shocks = zeros(nshocks, forecast_horizons)

    for t in 1:forecast_horizons
        shocks[:, t] = rand(dist)
    end

    compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
        shocks, z)
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
