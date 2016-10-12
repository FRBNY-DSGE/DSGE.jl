"""
```
smooth_all(m, data, systems, kals; cond_type = :none, procs = [myid()])
```

Computes and returns the smoothed values of states for every parameter draw.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `Matrix` of data for observables. This should
  include the conditional period if `cond_type={:semi,:full}`
- `systems::DVector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `kals::DVector{Kalman}`: a vector of `Kalman` objects containing the results of
  the Kalman filter for each draw

### Keyword Arguments

- `cond_type`: conditional case. See `forecast_all` for documentation of all cond_type options.`
- `procs`: list of worker processes over which to distribute draws. Defaults to `[myid()]`. Note that `

### Outputs

- `smoothed_states`: an `ndraws` x `nstates` x `nperiods` `DArray` of
  smoothed states returned from the smoother specified by
  `forecast_smoother(m)`.
- `smoothed_shocks`: an `ndraws` x `nshocks` x `nperiods` `DArray` of
  smoothed shocks returned from the smoother.
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
smooth(m, data, system, kals; cond_type = :none)
```

Computes and returns the smoothed values of states for the system `system`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `Matrix` of data for observables. This should
  include the conditional period if `cond_type={:semi,:full}`
- `system::System`: a `System` object representing the state-space system
- `kal::Kalman`: a `Kalman` object containing the results of
  the Kalman filter for this system

### Keyword Arguments

- `cond_type`: conditional case. See `forecast_all` for documentation of all cond_type options.`

### Outputs

- `smoothed_states`: an `ndraws` x `nstates` x `nperiods` `DArray` of
  smoothed states returned from the smoother specified by
  `forecast_smoother(m)`.
- `smoothed_shocks`: an `ndraws` x `nshocks` x `nperiods` `DArray` of
  smoothed shocks returned from the smoother.
"""
function smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System{S}, kal::Kalman{S})

    alpha_hat, eta_hat = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, system, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, system, kal[:z0], kal[:vz0])
    end

    return alpha_hat, eta_hat
end

function smooth{S<:AbstractFloat}(m::AbstractModel, data::DataFrame,
    system::System{S}, kal::Kalman{S})

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth(m, data, system, kal)
 end