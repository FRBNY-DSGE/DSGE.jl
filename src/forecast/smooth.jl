"""
```
smooth(m, df, system, kal; cond_type = :none, draw_states = true,
    include_presample = false)
```

Computes and returns the smoothed values of states and shocks for the system
`system`.

### Inputs

- `m::AbstractModel`: model object
- `df`: `DataFrame` of data for observables. This should include the conditional
  period if `cond_type in [:semi, :full]`
- `system::System`: `System` object representing the state-space system
- `kal::Kalman`: `Kalman` object containing the results of the Kalman filter for
  this system

### Keyword Arguments

- `cond_type`: conditional case. See `forecast_one` for documentation of all
  `cond_type` options
- `draw_states`: if using a simulation smoother (i.e.
  `forecast_smoother(m) in [:carter_kohn, :durbin_koopman]`), indicates whether
   to draw smoothed states from the distribution `N(z_{t|T}, P_{t|T})` or to use
   the mean `z_{t|T}`. Defaults to `false`. If not using a simulation smoother,
   this flag has no effect (though the user will be warned if
   `draw_states = true`)
- `include_presample::Bool`: indicates whether to include presample periods in
  the returned matrices. Defaults to `false`.
- `in_sample::Bool`: indicates whether or not to discard out of sample rows in `df_to_matrix` call.

### Outputs

- `states::Matrix{S}`: array of size `nstates` x `hist_periods` of smoothed
  states (not including the presample)
- `shocks::Matrix{S}`: array of size `nshocks` x `hist_nperiods` of smoothed
  shocks
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `hist_periods` of
  pseudo-observables computed from the smoothed states
- `initial_states::Vector{S}`: vector of length `nstates` of the smoothed states
  in the last presample period. This is used as the initial state for computing
  the deterministic trend

### Notes

`states` and `shocks` are returned from the smoother specified by
`forecast_smoother(m)`, which defaults to `:durbin_koopman`. This can be
overridden by calling

```
m <= Setting(:forecast_smoother, :koopman_smoother))
```

before calling `smooth`.
"""
function smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System{S}, kal::Kalman{S}; cond_type::Symbol = :none,
    draw_states::Bool = false, include_presample::Bool = false,
    in_sample::Bool = true)

    data = df_to_matrix(m, df; cond_type = cond_type, in_sample = in_sample)

    # Partition sample into pre- and post-ZLB regimes
    # Note that the post-ZLB regime may be empty if we do not impose the ZLB
    start_date = max(date_presample_start(m), df[1, :date])
    regime_inds = zlb_regime_indices(m, data, start_date)

    # Get system matrices for each regime
    TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs = zlb_regime_matrices(m, system, start_date)

    # Call smoother
    smoother = eval(Symbol(forecast_smoother(m), "_smoother"))

    states, shocks = if smoother == hamilton_smoother
        draw_states ? warn("$smoother called with draw_states = true") : nothing
        smoother(regime_inds, data, TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs,
            kal[:z0], kal[:vz0])
    elseif smoother == koopman_smoother
        draw_states ? warn("$smoother called with draw_states = true") : nothing
        smoother(regime_inds, data, TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs,
            kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif smoother in [carter_kohn_smoother, durbin_koopman_smoother]
        smoother(regime_inds, data, TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs,
            kal[:z0], kal[:vz0]; draw_states = draw_states)
    else
        error("Invalid smoother: $(forecast_smoother(m))")
    end

    # Index out last presample period, used to compute the deterministic trend
    t0 = n_presample_periods(m)
    t1 = index_mainsample_start(m)
    if include_presample
        initial_states = kal[:z0]
    else
        initial_states = states[:, t0]
        states = states[:, t1:end]
        shocks = shocks[:, t1:end]
    end

    # Map smoothed states to pseudo-observables
    pseudo = system[:ZZ_pseudo] * states .+ system[:DD_pseudo]

    return states, shocks, pseudo, initial_states
end
