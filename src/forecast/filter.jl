"""
```
filter_all(m::AbstractModel, data, systems, z0, vz0; cond_type = :none,
      lead = 0, allout = false, include_presample = true, procs = [myid()])
```

Distributes input arguments over `procs`; then computes and returns
the filtered values of states for every state-space system in
`systems`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `systems::DVector{System{S}}` vector of `System` objects specifying
  state-space system matrices for each draw
- `z0::Matrix{S}`: optional `Nz` x 1 initial state vector
- `vz0::Matrix{S}`: optional `Nz` x `Nz` initial state covariance matrix

where `S<:AbstractFloat`.

### Keyword Arguments

- `cond_type::Symbol`: conditional case. See `forecast_all` for documentation of
  all cond_type options.
- `lead::Int`: number of steps to forecast after the end of the data. Defaults
  to 0.
- `allout::Bool`: optional keyword argument indicating whether we want optional
  output variables returned as well. Defaults to `false`.
- `include_presample::Bool`: indicates whether to include presample periods in
  the returned vector of `Kalman` objects. Defaults to `true`.
- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`

### Outputs

- `kals::DVector{Kalman{S}`: vector of `ndraws` many `Kalman` objects. See
  `Kalman` documentation for more details.
"""
function filter_all{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    systems::DVector{System{S}, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
    include_presample::Bool = true, procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filter_all(m, data, systems, z0, vz0; lead = lead, allout = allout,
           include_presample = include_presample, procs = procs)
end

function filter_all{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    systems::DVector{System{S}, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0, allout::Bool = false, include_presample::Bool = true,
    procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)

    # Numbers of useful things
    ndraws = length(systems)
    nprocs = length(procs)

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), procs, [nprocs])
    datas  = dfill(data, (ndraws,), procs, [nprocs])
    z0s    = dfill(z0,   (ndraws,), procs, [nprocs])
    vz0s   = dfill(vz0,  (ndraws,), procs, [nprocs])

    # Construct distributed array of Kalman objects
    kals = DArray((ndraws,), procs, [nprocs]) do I
        draw_inds = first(I)
        ndraws_local = length(draw_inds)
        localpart = Vector{Kalman{S}}(ndraws_local)

        for i in draw_inds
            i_local = mod(i-1, ndraws_local) + 1
            localpart[i_local] = filter(models[i], datas[i], systems[i], z0s[i], vz0s[i];
                                        allout = allout, include_presample = include_presample)
        end
        return localpart
    end
end


"""
```
filter(m, data, system, z0, vz0; cond_type = :none, lead = 0, allout = false,
      include_presample = true)
```

Computes and returns the filtered values of states for the state-space
system corresponding to the current parameter values of model `m`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `system::System`: `System` object specifying state-space system matrices for
  the model
- `z0::Vector{S}`: optional `Nz` x 1 initial state vector
- `vz0::Matrix{S}`: optional `Nz` x `Nz` initial state covariance matrix

where `S<:AbstractFloat`.

### Keyword Arguments

- `cond_type::Symbol`: conditional case. See `forecast_all` for documentation of
  all `cond_type` options.
- `lead::Int`: number of steps to forecast after the end of the data. Defaults
  to 0.
- `allout::Bool`: optional keyword argument indicating whether we want optional
  output variables returned as well. Defaults to `false`.
- `include_presample::Bool`: indicates whether to include presample periods in
  the returned vector of `Kalman` objects. Defaults to `true`.

### Outputs

- `kal::Kalman`: see `Kalman` documentation for more details.
"""
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, system::System{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
    include_presample::Bool = true)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, system, z0, vz0; lead = lead, allout = allout,
           include_presample = include_presample)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, system::System{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0, allout::Bool = false, include_presample::Bool = true)

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]
    V_all   = system[:VVall]

    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, T, R, C, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    else
        # Regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, T, C, Z, D, V_all, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    end

    return kal
end

"""
```
filterandsmooth_all(m, data, systems, z0 = Vector{S}(), vz0 = Matrix{S}();
    cond_type = :none, lead = 0, procs = [myid()])
```

Computes and returns the smoothed states, shocks, and pseudo-observables, as
well as the filtered states for the last historical period, for every
state-space system in `systems`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `systems::DVector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw
- `z0::Vector{S}`: optional `Nz` x 1 initial state vector
- `vz0::Matrix{S}`: optional `Nz` x `Nz` initial state covariance matrix

where `S<:AbstractFloat`.

### Outputs

- `states::DArray{S, 3}`: array of size `ndraws` x `nstates` x `hist_periods` of
  smoothed states for each draw
- `shocks::DArray{S, 3}`: array of size `ndraws` x `nshocks` x `hist_periods` of
  smoothed shocks for each draw
- `pseudo::DArray{S, 3}`: array of size `ndraws` x `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
- `zends::DVector{Vector{S}}`: vector of `ndraws` many final state vectors (each
  length `nstates`) for the corresponding draws

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `filterandsmooth_all`.
"""
function filterandsmooth_all{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    systems::DVector{System{S}, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, lead::Int = 0,
    procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filterandsmooth_all(m, data, systems, z0, vz0; lead = lead, procs = procs)
end

function filterandsmooth_all{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    systems::DVector{System{S}, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0, procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)

    # Numbers of useful things
    ndraws = length(systems)
    nprocs = length(procs)
    nperiods = size(data, 2) - n_presample_periods(m)

    nstates = n_states_augmented(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    shocks_range = (nstates + 1):(nstates + nshocks)
    pseudo_range = (nstates + nshocks + 1):(nstates + nshocks + npseudo)
    zend_range   = nstates + nshocks + npseudo + 1

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), procs, [nprocs])
    datas  = dfill(data, (ndraws,), procs, [nprocs])
    z0s    = dfill(z0,   (ndraws,), procs, [nprocs])
    vz0s   = dfill(vz0,  (ndraws,), procs, [nprocs])

    # Construct distributed array of smoothed states, shocks, and pseudo-observables
    out = DArray((ndraws, nstates + nshocks + npseudo + 1, nperiods), procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = length(draw_inds)

        for i in draw_inds
            states, shocks, pseudo, zend = filterandsmooth(models[i], datas[i], systems[i], z0s[i], vz0s[i])

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = states
            localpart[i_local, shocks_range, :] = shocks
            if forecast_pseudoobservables(m)
                localpart[i_local, pseudo_range, :] = pseudo
            end
            localpart[i_local, zend_range,   1:nstates] = zend
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:nperiods])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:nperiods])
    zend   = DArray((ndraws,), procs, [nprocs]) do I
        Vector{S}[convert(Array, slice(out, i, zend_range, 1:nstates)) for i in first(I)]
    end

    # Index out SubArray for each smoothed type
    return states, shocks, pseudo, zend
end

"""
```
filterandsmooth(m, data, system, z0 = Vector{S}(), vz0 = Matrix{S}(); lead = 0)
```

Computes and returns the smoothed states, shocks, and pseudo-observables, as
well as the filtered states for the last historical period, for the state-space
system given by `system`.

### Inputs

- `m::AbstractModel`: model object
- `data`: `DataFrame` or `nobs` x `hist_periods` `Matrix{S}` of data for
  observables. This should include the conditional period if `cond_type in
  [:semi, :full]`
- `system::System{S}`: `System` objects specifying state-space system matrices
- `z0::Vector{S}`: optional `Nz` x 1 initial state vector
- `vz0::Matrix{S}`: optional `Nz` x `Nz` initial state covariance matrix

where `S<:AbstractFloat`.

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `hist_periods` of smoothed
  states
- `shocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of smoothed
  shocks
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
- `zend::Vector{S}`: final state vector of length `nstates`

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `filterandsmooth`.
"""
function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System{S}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0)

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]
    V_all   = system[:VVall]

    Z_pseudo, D_pseudo = if forecast_pseudoobservables(m)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    filterandsmooth(m, data, T, R, C, Q, Z, D, V_all, Z_pseudo, D_pseudo, z0, vz0; lead = lead)
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, V_all::Matrix{S},
    Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0)

    ## 1. Filter

    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, T, R, C, z0, vz0; lead =
            lead, allout = true, include_presample = true)
    else
        # Regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, T, C, Z, D, V_all, z0, vz0;
            lead = lead, allout = true, include_presample = true)
    end

    ## 2. Smooth

    states, shocks = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, T, R, C, Q, Z, D,
                        kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, T, R, C, Q, Z, D,
                                kal[:z0], kal[:vz0])
    end

    ## 3. Map smoothed states to pseudo-observables

    pseudo = if forecast_pseudoobservables(m)
        D_pseudo .+ Z_pseudo * states
    else
        Matrix{S}()
    end

    return states, shocks, pseudo, kal[:zend]
end
