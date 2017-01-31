"""
```
filter(m, data, system, z0, vz0; cond_type = :none, allout = false,
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
- `allout::Bool`: optional keyword argument indicating whether we want optional
  output variables returned as well. Defaults to `false`.
- `include_presample::Bool`: indicates whether to include presample periods in
  the returned vector of `Kalman` objects. Defaults to `true`.

### Outputs

- `kal::Kalman`: see `Kalman` documentation for more details.
"""
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, system::System{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, allout::Bool = false,
    include_presample::Bool = true)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, system, z0, vz0; allout = allout,
           include_presample = include_presample)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, system::System{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    allout::Bool = false, include_presample::Bool = true)

    # Unpack system
    T, R, C     = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D     = system[:QQ], system[:ZZ], system[:DD]
    M, E, V_all = system[:MM], system[:EE], system[:VVall]

    if n_anticipated_shocks(m) > 0
        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, T, R, C, z0, vz0;
                           ZZ = Z, DD = D, QQ = Q, MM = M, EE = E, VVall = V_all,
                           allout = true, include_presample = true)
    else
        # Regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, T, C, Z, D, V_all, z0, vz0;
            allout = allout, include_presample = include_presample)
    end

    return kal
end

"""
```
filterandsmooth(m, data, system, z0 = Vector{S}(), vz0 = Matrix{S}();
    cond_type = :none)
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

### Keyword Arguments

- `cond_type::Symbol`: conditional case. See `forecast_all` for documentation of
  all cond_type options.

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `hist_periods` of smoothed
  states
- `shocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of smoothed
  shocks
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :kalman_smoother))
```

before calling `filterandsmooth`.
"""
function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System{S}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filterandsmooth(m, data, system, z0, vz0)
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System{S}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}())

    # 0. Unpack system
    T, R, C     = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D     = system[:QQ], system[:ZZ], system[:DD]
    M, E, V_all = system[:MM], system[:EE], system[:VVall]

    Z_pseudo, D_pseudo = if forecast_pseudoobservables(m)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    ## 1. Filter

    if n_anticipated_shocks(m) > 0
        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, T, R, C, z0, vz0;
                           ZZ = Z, DD = D, QQ = Q, MM = M, EE = E, VVall = V_all,
                           allout = true, include_presample = true)
    else
        # Regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, T, C, Z, D, V_all, z0, vz0;
                  allout = true, include_presample = true)
    end

    ## 2. Smooth

    states, shocks = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, T, R, C, Q, Z, D,
                        kal[:z0], kal[:vz0], kal[:pred], kal[:vpred];
                        include_presample = true)
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, M, E, V_all,
                                kal[:z0], kal[:vz0];
                                include_presample = true)
    end

    # Index out last presample states, used to compute the deterministic trend
    t0 = n_presample_periods(m)
    t1 = index_mainsample_start(m)
    initial_states = states[:, t0]
    states = states[:, t1:end]
    shocks = shocks[:, t1:end]

    ## 3. Map smoothed states to pseudo-observables

    pseudo = if forecast_pseudoobservables(m)
        D_pseudo .+ Z_pseudo * states
    else
        Matrix{S}()
    end

    return states, shocks, pseudo, initial_states
end
