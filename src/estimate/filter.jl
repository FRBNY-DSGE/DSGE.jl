"""
```
filter(m, data, system, z0, P0; cond_type = :none, allout = true,
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
  output variables (besides `log_likelihood`, `zend`, and `Pend`) returned as
  well. Defaults to `true`.
- `include_presample::Bool`: indicates whether to include presample periods in
  the returned vector of `Kalman` objects. Defaults to `true`.

### Outputs

- `kal::Kalman`: see `Kalman` documentation for more details.
"""
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, system::System{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, allout::Bool = true,
    include_presample::Bool = true)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, system, z0, P0; allout = allout,
        include_presample = include_presample)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    catch_errors::Bool = false, allout::Bool = true,
    include_presample::Bool = true)

    # If we are in Metropolis-Hastings, then any errors coming out of `gensys`
    # should be caught and a -Inf posterior should be returned.
    system = try
        Nullable{System}(compute_system(m))
    catch err
        if catch_errors && isa(err, GensysError)
            warn(err.msg)
            Nullable{System}()
        else
            rethrow(err)
        end
    end

    if isnull(system)
        return Kalman(-Inf)
    else
        try
            filter(m, data, get(system), z0, P0; allout = allout,
                   include_presample = include_presample)
        catch err
            if catch_errors && isa(err, DomainError)
                warn("Log of incremental likelihood is negative; returning -Inf")
                return Kalman(-Inf)
            else
                rethrow(err)
            end
        end
    end
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, system::System,
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = true, include_presample::Bool = true)

    # Partition sample into pre- and post-ZLB regimes
    # Note that the post-ZLB regime may be empty if we do not impose the ZLB
    regime_inds = zlb_regime_indices(m, data)

    # Get system matrices for each regime
    TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs = zlb_regime_matrices(m, system)

    # If z0 and P0 provided, check that rows and columns corresponding to
    # anticipated shocks are zero in P0
    if !isempty(z0) && !isempty(P0)
        ant_state_inds = setdiff(1:n_states_augmented(m), inds_states_no_ant(m))
        @assert all(x -> x == 0, P0[:, ant_state_inds])
        @assert all(x -> x == 0, P0[ant_state_inds, :])
    end

    # Specify number of presample periods if we don't want to include them in
    # the final results
    T0 = include_presample ? 0 : n_presample_periods(m)

    # Run Kalman filter, construct Kalman object, and return
    out = kalman_filter(regime_inds, data, TTTs, RRRs, CCCs,
              QQs, ZZs, DDs, EEs, z0, P0;
              allout = allout, n_presample_periods = T0)

    return Kalman(out...)
end