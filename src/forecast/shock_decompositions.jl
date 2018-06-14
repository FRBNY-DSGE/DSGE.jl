"""
```
shock_decompositions(m, system, histshocks)

shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)
```

### Inputs

- `system::System{S}`: state-space system matrices
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
  pseudo-observable shock decompositions

where `nperiods = `end_index - start_index + 1`.
"""
function DSGE.shock_decompositions{S<:AbstractFloat}(m::AbstractModel,
    system::System{S}, histshocks::Matrix{S};
    shock_start_date::Date = date_mainsample_start(m),
    shock_end_date::Date = DSGE.iterate_quarters(date_presample_end(m), size(histshocks, 2)))

    horizon         = forecast_horizons(m)

    # Indices for which to return shockdecs
    keep_start_ind  = DSGE.index_shockdec_start(m)
    keep_end_ind    = DSGE.index_shockdec_end(m)

    # Indices for which to apply shocks while computing shockdecs
    shock_start_ind = DSGE.subtract_quarters(shock_start_date, date_presample_end(m))
    shock_end_ind   = DSGE.subtract_quarters(shock_end_date, date_presample_end(m))

    shock_decompositions(system, horizon, histshocks, keep_start_ind, keep_end_ind,
                         shock_start_index = shock_start_ind, shock_end_index = shock_end_ind)
end

function DSGE.shock_decompositions{S<:AbstractFloat}(system::System{S},
    forecast_horizons::Int, histshocks::Matrix{S},
    keep_start_index::Int, keep_end_index::Int;
    shock_start_index::Int = 1, shock_end_index = size(histshocks, 2))

    # Setup
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 2)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)
    histperiods  = size(histshocks, 2)
    allperiods   = histperiods + forecast_horizons

    states = zeros(S, nstates, allperiods, nshocks)
    obs    = zeros(S, nobs, allperiods, nshocks)
    pseudo = zeros(S, npseudo, allperiods, nshocks)

    # Check dates
    if forecast_horizons <= 0 || keep_start_index < 1 || keep_end_index > allperiods
        throw(DomainError())
    end

    # Set constant system matrices to 0
    system = DSGE.zero_system_constants(system)

    s_0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, allperiods)
        shocks[i, shock_start_index:shock_end_index] = histshocks[i, shock_start_index:shock_end_index]

        # Iterate state space forward
        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(system, s_0, shocks)
    end

    # Return shock decompositions in appropriate range
    if keep_start_index == 1 && keep_end_index == allperiods
        return states, obs, pseudo
    else
        range = keep_start_index:keep_end_index
        return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
    end
end

"""
```
deterministic_trends(m, system, s_0)

deterministic_trends(system, s_0, nperiods, start_index, end_index)
```

Compute deterministic trend values of states, observables, and
pseudo-observables, given a model object and system matrices. The deterministic
trend for a single draw is simply the series that would be obtained by iterating
the state-space system forward, beginning from a state vector `s_0` in the last
presample period.

### Inputs

- `m::AbstractModel`: model object
- `system::System{S}`: state-space system matrices
- `s_0`::Vector{S}: initial state vector

where `S<:AbstractFloat`.

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `nperiods` of state
  steady-state values
- `obs::Matrix{S}`: matrix of size `nobs` x `nperiods` of observable
  steady-state values
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `nperiods` of
  pseudo-observable steady-state values

where `nperiods` is the number of quarters between
`date_shockdec_start(m)` and `date_shockdec_end(m)`, inclusive.
"""
function DSGE.deterministic_trends{S<:AbstractFloat}(m::AbstractModel{S},
    system::System{S}, s_0::Vector{S})

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods    = subtract_quarters(date_forecast_end(m), date_mainsample_start(m)) + 1
    start_index = DSGE.index_shockdec_start(m)
    end_index   = DSGE.index_shockdec_end(m)

    deterministic_trends(system, s_0, nperiods, start_index, end_index)
end

function DSGE.deterministic_trends{S<:AbstractFloat}(system::System{S}, s_0::Vector{S}, nperiods::Int,
    start_index::Int, end_index::Int)

    # Set constant system matrices to 0
    system = DSGE.zero_system_constants(system)

    # Construct matrix of 0 shocks for entire history and forecast horizon
    nshocks  = size(system[:RRR], 2)
    shocks   = zeros(S, nshocks, nperiods)

    # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
    states, obs, pseudo, _ = forecast(system, s_0, shocks)

    # Return deterministic trends in appropriate range
    if start_index == 1 && end_index == nperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range], obs[:, range], pseudo[:, range]
    end
end

"""
```
trends{S<:AbstractFloat}(system::System{S})
```

Compute trend (steady-state) states, observables, and pseudo-observables. The
trend is used for plotting shock decompositions.
"""
function DSGE.trends{S<:AbstractFloat}(system::System{S})

    state_trend  = system[:CCC]
    obs_trend    = system[:ZZ]*system[:CCC] + system[:DD]
    pseudo_trend = system[:ZZ_pseudo]*system[:CCC] + system[:DD_pseudo]

    return state_trend, obs_trend, pseudo_trend
end