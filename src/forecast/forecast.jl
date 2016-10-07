"""
```
forecast{T<:AbstractFloat}(m::AbstractModel, system::Vector{System{T}},
    initial_state_draws::Vector{Vector{T}}; shock_distributions::Union{Distribution,
    Matrix{T}} = Matrix{T}())
```

Computes forecasts for all draws, given a model object, system matrices, and a
matrix of shocks or a distribution of shocks

### Inputs

- `m`: model object
- `systems::Vector{System}`: a vector of `System` objects specifying state-space
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
function forecast{S<:AbstractFloat}(m::AbstractModel,
    systems::DVector{System{S}, Vector{System{S}}},
    z0s::DVector{Vector{S}, Vector{Vector{S}}};
    shock_distributions::Union{Distribution,Matrix{S}} = Matrix{S}(),
    procs::Vector{Int} = [myid()])

    # Numbers of useful things
    ndraws = length(systems)
    nprocs = length(procs)
    horizon = forecast_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)
    shocks_range = (nstates + nobs + npseudo + 1):(nstates + nobs + npseudo + nshocks)

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
            dfill(zeros(nshocks, horizon), (ndraws,), procs, [nprocs])
        else
            # use t-distributed shocks
            if forecast_tdist_shocks(m)
                dfill(Distributions.TDist(forecast_tdist_df_val(m)), (ndraws,), procs, [nprocs])
            # use normally distributed shocks
            else
                DArray(I -> [DegenerateMvNormal(zeros(nshocks), sqrt(s[:QQ])) for s in systems[I...]],
                       (ndraws,), procs, [nprocs])
            end
        end
    end

    # Construct distributed array of forecast outputs
    out = DArray((ndraws, nstates + nobs + npseudo + nshocks, horizon), procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds
            # Get pseudomeasurement matrices
            Z_pseudo, D_pseudo = if forecast_pseudoobservables(m)
                _, pseudo_mapping = pseudo_measurement(m)
                pseudo_mapping.ZZ, pseudo_mapping.DD
            else
                Matrix{S}(), Vector{S}()
            end

            states, obs, pseudo, shocks = compute_forecast(systems[i], horizon,
                shock_distributions[i], z0s[i], Z_pseudo, D_pseudo)

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = states
            localpart[i_local, obs_range,    :] = obs
            localpart[i_local, pseudo_range, :] = pseudo
            localpart[i_local, shocks_range, :] = shocks
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:horizon])
    obs    = convert(DArray, out[1:ndraws, obs_range,    1:horizon])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:horizon])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:horizon])

    return states, obs, pseudo, shocks
end

"""
```
compute_forecast(T, R, C, Z, D, forecast_horizons,
    shocks, z, Z_pseudo, D_pseudo)
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

- `:states`
- `:observables`
- `:pseudo_observables`
- `:shocks`
"""
function compute_forecast{S<:AbstractFloat}(system::System{S},
    forecast_horizons::Int, shocks::Matrix{S}, z0::Vector{S},
    Z_pseudo::Matrix{S} = Matrix{S}(), D_pseudo::Vector{S} = Vector{S}())

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Z, D = system[:ZZ], system[:DD]

    compute_forecast(T, R, C, Z, D, forecast_horizons, shocks, z0, Z_pseudo, D_pseudo)
end

function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Z::Matrix{S}, D::Vector{S}, forecast_horizons::Int, shocks::Matrix{S},
    z0::Vector{S}, Z_pseudo::Matrix{S} = Matrix{S}(), D_pseudo::Vector{S} = Vector{S}())

    if forecast_horizons <= 0
        throw(DomainError())
    end

    # Setup
    nshocks = size(R, 2)
    nstates = size(T, 2)
    nobs    = size(Z, 1)
    npseudo = size(Z_pseudo, 1)

    states = zeros(nstates, forecast_horizons)

    # Define our iteration function
    iterate(z_t1, ϵ_t) = C + T*z_t1 + R*ϵ_t

    # Iterate state space forward
    states[:, 1] = iterate(z0, shocks[:, 1])
    for t in 2:forecast_horizons
        states[:, t] = iterate(states[:, t-1], shocks[:, t])
    end

    # Apply measurement and pseudo-measurement equations
    obs    = D .+ Z*states
    pseudo = if !isempty(Z_pseudo) && !isempty(D_pseudo)
        D_pseudo .+ Z_pseudo * states
    else
        Matrix{S}()
    end

    # Return forecasts
    return states, obs, pseudo, shocks
end

# Utility method to actually draw shocks
function compute_forecast{S<:AbstractFloat}(system::System{S},
    forecast_horizons::Int, dist::Distribution, z0::Vector{S},
    Z_pseudo::Matrix{S} = Matrix{S}(), D_pseudo::Vector{S} = Vector{S}())

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Z, D    = system[:ZZ], system[:DD]

    compute_forecast(T, R, C, Z, D, forecast_horizons, dist, z0, Z_pseudo, D_pseudo)
end

function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Z::Matrix{S}, D::Vector{S}, forecast_horizons::Int, dist::Distribution,
    z0::Vector{S}, Z_pseudo::Matrix{S} = Matrix{S}(), D_pseudo::Vector{S} = Vector{S}())

    if forecast_horizons <= 0
        throw(DomainError())
    end

    nshocks = size(R, 2)
    shocks = zeros(nshocks, forecast_horizons)

    for t in 1:forecast_horizons
        shocks[:, t] = rand(dist)
    end

    compute_forecast(T, R, C, Z, D, forecast_horizons, shocks, z0,
                     Z_pseudo, D_pseudo)
end
