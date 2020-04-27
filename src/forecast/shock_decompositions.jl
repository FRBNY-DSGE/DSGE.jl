"""
```
shock_decompositions(m, system, histshocks)

shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)

shock_decompositions(m, system, histshocks, start_date, end_date)

shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index, regime_inds)
```

### Inputs

- `system::System{S}` or `RegimeSwitchingSystem{S}`: state-space system matrices
- `forecast_horizons::Int`: number of periods ahead to forecast
- `histshocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of
  historical smoothed shocks
- `start_index::Int`: first index from which to return computed shock
  decompositions
- `end_index::Int`: last index for which to return computed shock decompositions
- `start_date::Date`: initial date for deterministic trends
- `end_date::Date`: final date for deterministic trends
- `regime_inds::Vector{UnitRange}`: indices of the data corresponding to each regime.

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
function shock_decompositions(m::AbstractDSGEModel,
    system::System{S}, histshocks::Matrix{S}) where {S<:AbstractFloat}

    horizon   = forecast_horizons(m)
    start_ind = index_shockdec_start(m)
    end_ind   = index_shockdec_end(m)

    shock_decompositions(system, horizon, histshocks, start_ind, end_ind)
end

function shock_decompositions(system::System{S},
    forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int) where {S<:AbstractFloat}

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
    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError())
    end

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    z0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, allperiods)
        shocks[i, 1:histperiods] = histshocks[i, :]

        # Iterate state space forward
        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(system, z0, shocks)
    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
    end
end

function shock_decompositions(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, histshocks::Matrix{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m))) where {S<:AbstractFloat}
    # Set up dates
    horizon     = forecast_horizons(m)
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Set up regime indices. Note that we do not need to account
    # for ZLB split b/c `histshocks` should have zeros for anticipated shocks
    # in the pre-ZLB periods.
    regime_inds = regime_indices(m, start_date, end_date)

    shock_decompositions(system, horizon, histshocks, start_index, end_index, regime_inds)
end

function shock_decompositions(system::RegimeSwitchingSystem{S},
                              forecast_horizons::Int, histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}}) where {S<:AbstractFloat}

    # Setup
    nshocks     = size(system[1, :RRR], 2)
    nstates     = size(system[1, :TTT], 2)
    nobs        = size(system[1, :ZZ], 1)
    npseudo     = size(system[1, :ZZ_pseudo], 1)
    histperiods = size(histshocks, 2)
    allperiods  = histperiods + forecast_horizons

    states = zeros(S, nstates, allperiods, nshocks)
    obs    = zeros(S, nobs,    allperiods, nshocks)
    pseudo = zeros(S, npseudo, allperiods, nshocks)
    shocks = zeros(S, nshocks, histperiods)

    # Check dates
    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods || regime_inds[1][1] != 1 || regime_inds[end][end] != histperiods
        throw(DomainError())
    end

    # Set constant system matrices to 0 and start at stationary steady state
    system = zero_system_constants(system)
    z0     = zeros(S, nstates)

    fcast_shocks = zeros(S, nshocks, forecast_horizons) # these are always zero
    for i = 1:nshocks
        # Isolate single shock
        shocks[i, :] = histshocks[i, :]

        # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
        for (reg_num, reg_ind) in enumerate(regime_inds)
            init_state = (reg_num == regimes[1]) ? zeros(S, nstates) : states[:, reg_ind[1] - 1, i] # update initial state

            # Construct matrix of 0 shocks for this regime
            shocks = zeros(S, nshocks, length(reg_ind))
            sys    = system[reg_num] # grab the correct system
            states[:, reg_ind, i], obs[:, reg_ind, i], pseudo[:, reg_ind, i], _ = forecast(system, init_state, shocks[:, reg_ind])
        end

        fcast_inds = regime_inds[end][end] + 1:allperiods
        states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
            forecast(system[length(regime_inds)], states[:, regime_inds[end][end], i], fcast_shocks)

        # zero out shocks
        shocks[i, 1:histperiods] .= 0.
    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
    end
end

"""
```
deterministic_trends(m, system, z0)

deterministic_trends(system, z0, nperiods, start_index, end_index)

deterministic_trends(m, system, z0, start_date, end_date)

deterministic_trends(system, z0, nperiods, regime_inds, regimes)
```

Compute deterministic trend values of states, observables, and
pseudo-observables, given a model object and system matrices. The deterministic
trend for a single draw is simply the series that would be obtained by iterating
the state-space system forward, beginning from a state vector `z0` in the last
presample period.

### Inputs

- `m::AbstractDSGEModel`: model object
- `system::System{S}` or `RegimeSwitchingSystem`: state-space system matrices
- `z0`::Vector{S}: initial state vector
- `start_date::Date`: initial date for deterministic trends
- `end_date::Date`: final date for deterministic trends
- `regime_inds::Vector{UnitRange}`: indices of the data corresponding to each regime.
- `regimes::UnitRange`: which regimes are involved in the date range for which
    we want to compute the deterministic trends

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
function deterministic_trends(m::AbstractDSGEModel{S},
                              system::System{S}, z0::Vector{S}) where {S<:AbstractFloat}

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods    = subtract_quarters(date_forecast_end(m), date_mainsample_start(m)) + 1
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    deterministic_trends(system, z0, nperiods, start_index, end_index)
end

function deterministic_trends(system::System{S}, z0::Vector{S}, nperiods::Int,
                              start_index::Int, end_index::Int) where {S<:AbstractFloat}

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    # Construct matrix of 0 shocks for entire history and forecast horizon
    nshocks  = size(system[:RRR], 2)
    shocks   = zeros(S, nshocks, nperiods)

    # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
    states, obs, pseudo, _ = forecast(system, z0, shocks)

    # Return deterministic trends in appropriate range
    if start_index == 1 && end_index == nperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range], obs[:, range], pseudo[:, range]
    end
end

function deterministic_trends(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, z0::Vector{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m))) where {S<:AbstractFloat}

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods    = subtract_quarters(date_forecast_end(m), date_mainsample_start(m)) + 1
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Now intersect these periods with regime indices
    start_date  = max(date_presample_start(m), df[1, :date])
    regime_inds = regime_indices(m, start_date, end_date) # do not need to account for ZLB split b/c shocks don't matter for trends

    first_regime = findfirst(map(x -> start_index in x, regime_inds))
    last_regime  = findfirst(map(x -> end_index in x,   regime_inds))
    regimes      = first_regime:last_regime                  # need to specify which regimes we need for shockdec.
    regime_inds  = regime_inds[regimes]                      # also need to specify how long in each regime we are.
    regime_inds[1]   = start_index:regime_inds[1][end]       # if the start index is in the middle of a regime.
    regime_inds[end] = end_index:regime_inds[end][end]       # if the end index is in the middle of a regime.

    deterministic_trends(system, z0, nperiods, regime_inds, regimes)
end


function deterministic_trends(system::RegimeSwitchingSystem{S}, z0::Vector{S}, nperiods::Int,
                              regime_inds::Vector{UnitRange{Int}}, regimes::UnitRange{Int}) where {S<:AbstractFloat}

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    # Initialize storage matrices
    nshocks = size(system[1, :RRR], 2) # this is used in the forecasting loops to create shock matrices
    nstates = size(system[1, :TTT], 1)
    nobs    = size(system[1, :ZZ], 1)
    npseudo = size(system[1, :ZZ_pseudo], 1)
    states = Matrix{S}(undef, nstates, nperiods)
    obs    = Matrix{S}(undef, nobs,    nperiods)
    pseudo = Matrix{S}(undef, npseudo, nperiods)

    # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
    for (reg_num, reg_ind) in zip(regimes, regime_inds)
        init_state = (reg_num == regimes[1]) ? z0 : states[:, reg_ind[1] - 1] # update initial state

        # Construct matrix of 0 shocks for this regime
        shocks = zeros(S, nshocks, length(reg_ind))
        sys = system[reg_num] # grab the correct system
        states[:, reg_ind], obs[:, reg_ind], pseudo[:, reg_ind], _ = forecast(system, init_state, shocks)
    end

    return state, obs, pseudo
end

"""
```
trends(system::System{S}) where {S<:AbstractFloat}
trends(system::RegimeSwitchingSystem{S}) where {S<:AbstractFloat}
```

Compute trend (steady-state) states, observables, and pseudo-observables. The
trend is used for plotting shock decompositions.
"""
function trends(system::System{S}) where {S<:AbstractFloat}

    state_trend  = system[:CCC]
    obs_trend    = system[:ZZ]*system[:CCC] + system[:DD]
    pseudo_trend = system[:ZZ_pseudo]*system[:CCC] + system[:DD_pseudo]

    return state_trend, obs_trend, pseudo_trend
end

function trends(system::RegimeSwitchingSystem{S}) where {S<:AbstractFloat}

    nreg  = n_regimes(system)
    input = 1:nreg

    state_trends  = Matrix{S}(undef, length(system[1, :CCC]), nreg)
    obs_trends    = Matrix{S}(undef, length(system[1, :DD]), nreg)
    pseudo_trends = Matrix{S}(undef, length(system[1, :DD_pseudo]), nreg)
    for i in input
        state_trends[:, i]  = system[i, :CCC]
        obs_trends[:, i]    = system[i, :ZZ] * system[i, :CCC] + system[i, :DD]
        pseudo_trends[:, i] = system[i, :ZZ_pseudo] * system[i, :CCC] + system[i, :DD_pseudo]
    end

    return state_trends, obs_trends, pseudo_trends
end
