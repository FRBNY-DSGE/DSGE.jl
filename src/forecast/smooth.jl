"""
```
smooth_all(m, data, systems, kals; cond_type = :none)
```

Computes and returns the smoothed values of states for every parameter draw.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `systems::Vector{System{S}}`: vector of `ndraws` many `System` objects
  specifying state-space system matrices for each draw
- `kals::Vector{Kalman{S}}`: a vector of `ndraws` many `Kalman` objects
  containing the results of the Kalman filter for each draw

### Keyword Arguments

- `cond_type::Symbol`: conditional case. See `forecast_all` for documentation of
  all `cond_type` options.

### Outputs

- `states::Array{S, 3}`: array of size `ndraws` x `nstates` x `hist_periods` of
  smoothed states for each draw
- `shocks::Array{S, 3}`: array of size `ndraws` x `nshocks` x `hist_periods` of
  smoothed shocks for each draw
- `pseudo::Array{S, 3}`: array of size `ndraws` x `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
- `initial_states::Vector{Vector{{S}}`: vector of length `ndraws`, where each
  element is a vector of smoothed states in the last presample period. These are
  used as the initial state for computing the deterministic trend.

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `smooth_all`.
"""
function smooth_all{S<:AbstractFloat}(m::AbstractModel{S}, df::DataFrame,
    systems::Vector{System{S}}, kals::Vector{Kalman{S}};
    cond_type::Symbol = :none)

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth_all(m, data, systems, kals)
end

function smooth_all{S<:AbstractFloat}(m::AbstractModel{S}, data::Matrix{S},
    systems::Vector{System{S}}, kals::Vector{Kalman{S}})

    ndraws = length(systems)
    @assert length(kals) == ndraws

    nstates = n_states_augmented(m)
    nshocks = n_shocks_exogenous(m)
    npseudo = n_pseudoobservables(m)
    nperiods = size(data, 2) - n_presample_periods(m)

    states = zeros(S, ndraws, nstates, nperiods)
    shocks = zeros(S, ndraws, nshocks, nperiods)
    pseudo = zeros(S, ndraws, npseudo, nperiods)
    initial_states = Vector{Vector{S}}(ndraws)

    for i = 1:ndraws
        states_i, shocks_i, pseudo_i, initial_states_i = smooth(m, data, systems[i], kals[i])

        states[i, :, :] = states_i
        shocks[i, :, :] = shocks_i
        pseudo[i, :, :] = pseudo_i
        initial_states[i] = initial_states_i
    end

    return states, shocks, pseudo, initial_states
end



"""
```
smooth(m, data, system, kal; cond_type = :none)
```

Computes and returns the smoothed values of states and shocks for the system
`system`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `system::System`: `System` object representing the state-space system
- `kal::Kalman`: `Kalman` object containing the results of the Kalman filter for
  this system

### Keyword Arguments

- `cond_type`: conditional case. See `forecast_all` for documentation of all
  `cond_type` options.

### Outputs

- `states::Matrix{S}`: array of size `nstates` x `hist_periods` of smoothed
  states (not including the presample)
- `shocks::Matrix{S}`: array of size `nshocks` x `hist_nperiods` of smoothed
  shocks
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
- `initial_states::Vector{S}`: vector of length `nstates` of the smoothed states
  in the last presample period. This is used as the initial state for computing
  the deterministic trend

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `smooth`.
"""
function smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System{S}, kal::Kalman{S}; cond_type::Symbol = :none)

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth(m, data, system, kal)
 end

function smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System{S}, kal::Kalman{S})

    states, shocks = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, system, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred];
            include_presample = true)
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, system, kal[:z0], kal[:vz0];
            include_presample = true)
    end

    # Index out last presample states, used to compute the deterministic trend
    t0 = n_presample_periods(m)
    t1 = index_mainsample_start(m)
    initial_states = states[:, t0]
    states = states[:, t1:end]
    shocks = shocks[:, t1:end]

    # Map smoothed states to pseudo-observables
    pseudo = if forecast_pseudoobservables(m)
        system[:DD_pseudo] .+ system[:ZZ_pseudo] * states
    else
        Matrix{S}()
    end

    return states, shocks, pseudo, initial_states
end