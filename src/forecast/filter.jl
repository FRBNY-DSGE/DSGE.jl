"""
```
filter(m, data, system, z0, P0; cond_type = :none, allout = false,
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
- `P0::Matrix{S}`: optional `Nz` x `Nz` initial state covariance matrix

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
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none)

    data = df_to_matrix(m, df; cond_type = cond_type)

    kalman_filter(m, data, system, z0, P0;
        allout = true, include_presample = true)
end