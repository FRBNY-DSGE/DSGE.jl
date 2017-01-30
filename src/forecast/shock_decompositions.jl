"""
```
shock_decompositions(m, systems, histshocks)
```

Computes shock decompositions for all draws, given a model object, system
matrices, and historical smoothed shocks.

### Inputs

- `m::AbstractModel`: model object
- `systems::Vector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw
- `histshocks`::Array{S, 3}`: array of size `ndraws` x `nshocks` x
  `hist_periods` of smoothed historical shocks for each draw

where `S<:AbstractFloat`.

### Outputs

- `states::Array{S, 4}`: array of size `ndraws` x `nstates` x `nperiods` x
  `nshocks` of state shock decompositions for each draw
- `obs::Array{S, 4}`: array of size `ndraws` x `nobs` x `nperiods` x `nshocks`
  of observable shock decompositions for each draw
- `pseudo::Array{S, 4}`: array of size `ndraws` x `npseudo` x `nperiods` x
  `nshocks` of pseudo-observable shock decompositions for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.

where `nperiods` is the number of quarters between `date_shockdec_start(m)` and
`date_shockdec_end(m)`, inclusive. If `date_shockdec_start(m)` is null, shock
decompositions are returned beginning from `date_mainsample_start(m)`. Likewise, if
`date_shockdec_end(m)` is null, shock decompositions are returned up to
`date_forecast_end(m)`.
"""
function shock_decompositions{S<:AbstractFloat}(m::AbstractModel{S},
    systems::Vector{System{S}}, histshocks::Array{S, 3})

    # Numbers of useful things
    ndraws = length(systems)
    horizon = forecast_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    # Determine periods for which to return shock decompositions
    start_ind = index_shockdec_start(m)
    end_ind   = index_shockdec_end(m)
    nperiods  = end_ind - start_ind + 1

    states = zeros(ndraws, nstates, nperiods, nshocks)
    obs    = zeros(ndraws, nobs,    nperiods, nshocks)
    pseudo = zeros(ndraws, npseudo, nperiods, nshocks)

    for i = 1:ndraws
        shocks_i = squeeze(histshocks[i, :, :], 1)
        states_i, obs_i, pseudo_i = compute_shock_decompositions(systems[i], horizon,
            shocks_i, start_ind, end_ind)

        states[i, :, :, :] = states_i
        obs[i,    :, :, :] = obs_i
        pseudo[i, :, :, :] = pseudo_i
    end

    return states, obs, pseudo
end

"""
```
compute_shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)

compute_shock_decompositions(T, R, Z, D, Z_pseudo, D_pseudo, z0, shocks,
    forecast_horizons, histshocks, start_index, end_index)
```

### Inputs

- `system::System{S}`: state-space system matrices. Alternatively, provide
  transition equation matrices `T`, `R`; measurement equation matrices `Z`, `D`;
  and (possibly empty) pseudo-measurement equation matrices `Z_pseudo` and
  `D_pseudo`.
- `forecast_horizons::Int`: number of periods ahead to forecast
- `histshocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of
  historical smoothed shocks
- `start_index::Int`: first index from which to return computed shock
  decompositions
- `end_index::Int`: last index for which to return computed shock decompositions

where `S<:AbstractFloat`.

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `nperiods` x `nshocks` of
  state shock decompositions
- `obs::Array{S, 3}`: matrix of size `nobs` x `nperiods` x `nshocks` of
  observable shock decompositions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `nperiods` x `nshocks` of
  pseudo-observable shock decompositions. If the provided `Z_pseudo` and
  `D_pseudo` matrices are empty, then `pseudo` will be empty.

where `nperiods = `end_index - start_index + 1`.
"""
function compute_shock_decompositions{S<:AbstractFloat}(system::System{S},
    forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int)

    # Unpack system
    T, R = system[:TTT], system[:RRR]
    Z, D = system[:ZZ], system[:DD]

    Z_pseudo, D_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    compute_shock_decompositions(T, R, Z, D, Z_pseudo, D_pseudo,
        forecast_horizons, histshocks, start_index, end_index)
end

function compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S},
    R::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int)

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobs         = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    histperiods  = size(histshocks, 2)
    allperiods   = histperiods + forecast_horizons

    forecast_pseudo = !isempty(Z_pseudo) && !isempty(D_pseudo)

    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError())
    end

    states = zeros(S, nstates,      allperiods, nshocks)
    obs    = zeros(S, nobs, allperiods, nshocks)
    pseudo = if forecast_pseudo
        zeros(S, npseudo, allperiods, nshocks)
    else
        zeros(S, 0, 0, 0)
    end

    # Define our iteration function
    iterate(z_t1, ϵ_t) = T*z_t1 + R*ϵ_t
    z0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, allperiods)
        shocks[i, 1:histperiods] = histshocks[i, :]

        # Iterate state space forward
        states[:, 1, i] = iterate(z0, shocks[:, 1])
        for t in 2:allperiods
            states[:, t, i] = iterate(states[:, t-1, i], shocks[:, t])
        end

        # Apply measurement and pseudo-measurement equations
        obs[:, :, i] = D .+ Z * states[:, :, i]
        if forecast_pseudo
            pseudo[:, :, i] = D_pseudo .+ Z_pseudo * states[:, :, i]
        end

    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        if forecast_pseudo
            return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
        else
            return states[:, range, :], obs[:, range, :], pseudo
        end
    end
end
