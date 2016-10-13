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
    shocks::DArray{S, 3} = dzeros(S, (0, 0, 0), [myid()]),
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

    shocks_provided = !isempty(shocks)

    # Construct distributed array of forecast outputs
    out = DArray((ndraws, nstates + nobs + npseudo + nshocks, horizon), procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds
            # Index out shocks for draw i
            shocks_i = if shocks_provided
                convert(Array, slice(shocks, i, :, :))
            else
                Matrix{S}()
            end

            states, obs, pseudo, shocks = compute_forecast(m, systems[i],
                z0s[i]; shocks = shocks_i)

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
function compute_forecast{S<:AbstractFloat}(m::AbstractModel, system::System{S},
    z0::Vector{S}; shocks::Matrix{S} = Matrix{S}())

    # Numbers of things
    nshocks = n_shocks_exogenous(m)
    horizon = forecast_horizons(m)

    # Populate shocks matrix
    if isempty(shocks)
        shocks = zeros(nshocks, horizon)

        # Draw shocks if necessary
        if !forecast_kill_shocks(m)
            dist = if forecast_tdist_shocks(m)
                # Use t-distributed shocks
                Distributions.TDist(forecast_tdist_df_val(m))
            else
                # Use normally distributed shocks
                DegenerateMvNormal(zeros(nshocks), sqrt(system[:QQ]))
            end

            for t in 1:horizon
                shocks[:, t] = rand(dist)
            end
        end
    end

    compute_forecast(system, z0, shocks)
end


function compute_forecast{S<:AbstractFloat}(system::System{S}, z0::Vector{S},
    shocks::Matrix{S})

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]

    Z_pseudo, D_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    compute_forecast(T, R, C, Q, Z, D, Z_pseudo, D_pseudo, z0, shocks)
end

function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S},
    C::Vector{S}, Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, z0::Vector{S}, shocks::Matrix{S})

    # Setup
    nshocks = size(R, 2)
    nstates = size(T, 2)
    nobs    = size(Z, 1)
    npseudo = size(Z_pseudo, 1)
    horizon = size(shocks, 2)

    # Define our iteration function
    iterate(z_t1, ϵ_t) = C + T*z_t1 + R*ϵ_t

    # Iterate state space forward
    states = zeros(nstates, horizon)
    states[:, 1] = iterate(z0, shocks[:, 1])
    for t in 2:horizon
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
