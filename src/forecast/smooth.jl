"""
```
smooth_all(m, data, systems, kals; cond_type = :none, procs = [myid()])
```

Computes and returns the smoothed values of states for every parameter draw.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `systems::DVector{System{S}}`: vector of `ndraws` many `System` objects
  specifying state-space system matrices for each draw
- `kals::DVector{Kalman{S}}`: a vector of `ndraws` many `Kalman` objects
  containing the results of the Kalman filter for each draw

### Keyword Arguments

- `cond_type::Symbol`: conditional case. See `forecast_all` for documentation of
  all `cond_type` options.
- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`.

### Outputs

- `states::DArray{S, 3}`: array of size `ndraws` x `nstates` x `hist_periods` of
  smoothed states for each draw
- `shocks::DArray{S, 3}`: array of size `ndraws` x `nshocks` x `hist_periods` of
  smoothed shocks for each draw

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `smooth_all`.
"""
function smooth_all{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    systems::DVector{System{S}, Vector{System{S}}},
    kals::DVector{Kalman{S}}; cond_type::Symbol = :none,
    procs::Vector{Int} = [myid()])

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth_all(m, data, systems, kals; procs = procs)
end

function smooth_all{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    systems::DVector{System{S}, Vector{System{S}}},
    kals::DVector{Kalman{S}}; procs::Vector{Int} = [myid()])

    # numbers of useful things
    nprocs = length(procs)
    ndraws = length(systems)
    @assert length(kals) == ndraws

    nstates = n_states_augmented(m)
    nshocks = n_shocks_exogenous(m)
    nperiods = size(data, 2) - n_presample_periods(m)

    states_range = 1:nstates
    shocks_range = (nstates + 1):(nstates + nshocks)

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), procs, [nprocs])
    datas  = dfill(data, (ndraws,), procs, [nprocs])

    # Construct distributed array of smoothed states and shocks
    out = DArray((ndraws, nstates + nshocks, nperiods), procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = length(draw_inds)

        for i in draw_inds
            states, shocks = smooth(models[i], datas[i], systems[i], kals[i])

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = states
            localpart[i_local, shocks_range, :] = shocks
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:nperiods])

    return states, shocks
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
  states
- `shocks::Matrix{S}`: array of size `nshocks` x `hist_nperiods` of smoothed
  shocks

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

    alpha_hat, eta_hat = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, system, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, system, kal[:z0], kal[:vz0])
    end

    return alpha_hat, eta_hat
end