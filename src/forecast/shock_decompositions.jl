"""
```
shock_decompositions(m, system, histshocks)

shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)

shock_decompositions(m, system, histshocks, start_date, end_date)

shock_decompositions(m, system, forecast_horizons, histshocks, start_index,
    end_index, regime_inds, cond_type)
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
        throw(DomainError("At least one of forecast_horizons, the index for the start date of the shock decomposition, and the index for the end date of the shock decomposition is outside the permissible range."))
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
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none) where {S<:AbstractFloat}
    # Set up dates
    horizon     = forecast_horizons(m; cond_type = cond_type)
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Set up regime indices. Note that we do not need to account
    # for ZLB split b/c `histshocks` should have zeros for anticipated shocks
    # in the pre-ZLB periods.
    regime_inds = regime_indices(m, start_date, end_date)
    if regime_inds[1][1] < 1 # remove periods occuring before desired start date
        regime_inds[1] = 1:regime_inds[1][end]
    end

    shock_decompositions(m, system, horizon, histshocks, start_index, end_index, regime_inds, cond_type)
end

function shock_decompositions(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                              forecast_horizons::Int, histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol) where {S<:AbstractFloat}

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
    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError("At least one of forecast_horizons, the index for the start date of the shock decomposition, and the index for the end date of the shock decomposition is outside the permissible range."))
    end

    # Set constant system matrices to 0 and start at stationary steady state
    system = zero_system_constants(system)
    z0     = zeros(S, nstates)

    fcast_shocks = zeros(S, nshocks, forecast_horizons) # these are always zero
    for i = 1:nshocks
        # Isolate single shock
        shocks[i, :] = histshocks[i, :]

        # Looping through each historical regime
        for (reg_num, reg_ind) in enumerate(regime_inds)
            if reg_ind[end] > histperiods
                break
            end
            init_state = (reg_num == 1) ? zeros(S, nstates) : states[:, reg_ind[1] - 1, i] # update initial state

            states[:, reg_ind, i], obs[:, reg_ind, i], pseudo[:, reg_ind, i], _ =
                forecast(system[reg_num], init_state, shocks[:, reg_ind])
        end

        # Run forecast using regime-switching forecast
        if allperiods > histperiods
            fcast_inds = (histperiods + 1):allperiods
            states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
                forecast(m, system, states[:, histperiods, i], fcast_shocks, cond_type = cond_type)
        end

        # zero out shocks b/c want effects of single shocks
        shocks[i, :] .= 0.
    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
    end
end

# Get Shock Decompositions with an old system and a bar that gives
## effet of using new system compared to old.
function shock_decompositions(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, old_system::Int,
                              histshocks::Matrix{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none) where {S<:AbstractFloat}
    # Set up dates
    horizon     = forecast_horizons(m; cond_type = cond_type)
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Set up regime indices. Note that we do not need to account
    # for ZLB split b/c `histshocks` should have zeros for anticipated shocks
    # in the pre-ZLB periods.
    regime_inds = regime_indices(m, start_date, end_date)
    if regime_inds[1][1] < 1 # remove periods occuring before desired start date
        regime_inds[1] = 1:regime_inds[1][end]
    end

    shock_decompositions(m, system, old_system, horizon, histshocks, start_index, end_index, regime_inds, cond_type)
end

function shock_decompositions(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                              old_system::Int,
                              forecast_horizons::Int, histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol) where {S<:AbstractFloat}

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

    # Old System Values
    old_states = zeros(S, nstates, allperiods, nshocks)
    old_obs    = zeros(S, nobs,    allperiods, nshocks)
    old_pseudo = zeros(S, npseudo, allperiods, nshocks)

    # Check dates
    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError("At least one of forecast_horizons, the index for the start date of the shock decomposition, and the index for the end date of the shock decomposition is outside the permissible range."))
    end

    # Set constant system matrices to 0 and start at stationary steady state
    system = zero_system_constants(system)
    z0     = zeros(S, nstates)

    fcast_shocks = zeros(S, nshocks, forecast_horizons) # these are always zero
    for i = 1:nshocks
        # Isolate single shock
        shocks[i, :] = histshocks[i, :]

        # Old System Values
        # init_state_old = zeros(S, nstates)
        # old_states[:,1:histperiods,i], old_obs[:,1:histperiods,i], old_pseudo[:,1:histperiods,i], _ =
        #     forecast(old_system, init_state_old, shocks[:, 1:histperiods])

        counts = 1
        # Looping through each historical regime
        for (reg_num, reg_ind) in enumerate(regime_inds)
            if reg_ind[end] > histperiods
                break
            end
            init_state = (reg_num == 1) ? zeros(S, nstates) : states[:, reg_ind[1] - 1, i] # update initial state

            states[:, reg_ind, i], obs[:, reg_ind, i], pseudo[:, reg_ind, i], _ =
                forecast(system[reg_num], init_state, shocks[:, reg_ind])

            # Old System Values
            if reg_ind[end] < old_system
                old_system2 = system[reg_num]
                counts = reg_num
            else
                old_system2 = system[1]#system[counts]
            end
            init_state_old = (reg_num == 1) ? zeros(S, nstates) : old_states[:, reg_ind[1] - 1, i]
            old_states[:, reg_ind, i], old_obs[:, reg_ind, i], old_pseudo[:, reg_ind, i], _ =
                forecast(old_system2, init_state_old, shocks[:,reg_ind])
        end

        # Run forecast using regime-switching forecast
        if allperiods > histperiods
            fcast_inds = (histperiods + 1):allperiods
            states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
                forecast(m, system, states[:, histperiods, i], fcast_shocks, cond_type = cond_type)

            # Old System Values
            old_system2 = system[1]#system[counts]
            old_states[:, fcast_inds, i], old_obs[:, fcast_inds, i], old_pseudo[:, fcast_inds, i], _ =
                forecast(old_system2, old_states[:, histperiods, i], fcast_shocks)
        end

        # zero out shocks b/c want effects of single shocks
        shocks[i, :] .= 0.
    end

    # Add new system as shock to old system
    old_states = cat(old_states, zeros(size(old_states[:,:,1])), dims = 3)

    for i in 1:size(states,1)
        for j in 1:size(states,2)
            old_states[i,j,end] = sum([states[i,j,k] for k in 1:size(states,3)]) - sum([old_states[i,j,k] for k in 1:size(states,3)])
        end
    end

    for i in 1:size(obs,1)
        for j in 1:size(obs,2)
            old_obs[i,j,end] = sum([obs[i,j,k] for k in 1:size(obs,3)]) - sum([old_obs[i,j,k] for k in 1:size(obs,3)])
        end
    end

    for i in 1:size(pseudo,1)
        for j in 1:size(pseudo,2)
            old_pseudo[i,j,end] = sum([pseudo[i,j,k] for k in 1:size(pseudo,3)]) - sum([old_pseudo[i,j,k] for k in 1:size(pseudo,3)])
        end
    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return old_states, old_obs, old_pseudo
    else
        range = start_index:end_index
        return old_states[:, range, :], old_obs[:, range, :], old_pseudo[:, range, :]
    end
end


"""
```
deterministic_trends(m, system, z0)

deterministic_trends(system, z0, nperiods, start_index, end_index)

deterministic_trends(m, system, z0, start_date, end_date)

deterministic_trends(m, system, z0, nperiods, start_index, end_index,
    regime_inds, regimes, cond_type)
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
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none) where {S<:AbstractFloat}

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods    = subtract_quarters(date_forecast_end(m), date_mainsample_start(m)) + 1
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Now intersect these periods with regime indices
    regime_inds = regime_indices(m, start_date, end_date) # do not need to account for ZLB split b/c shocks don't matter for trends
    if regime_inds[1][1] < 1
        regime_inds[1] = 1:regime_inds[1][end]
    end
    if regime_inds[end][end] > end_index
        regime_inds[end] = regime_inds[end][1]:end_index # if the end index is in the middle of a regime or is past the regime's end
    end
    if length(regime_inds[end]) == 0
        pop!(regime_inds)
    end

    return deterministic_trends(m, system, z0, nperiods, start_index, end_index, regime_inds, cond_type)
end


function deterministic_trends(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S}, z0::Vector{S}, nperiods::Int,
                              start_index::Int, end_index::Int, regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol) where {S<:AbstractFloat}

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

    # Use forecast to iterate state-space system forward
    # without shocks or the ZLB procedure during history
    for (reg_num, reg_ind) in enumerate(regime_inds)
        init_state = (reg_num == 1) ? z0 : states[:, reg_ind[1] - 1] # update initial state

        # Construct matrix of 0 shocks for this regime
        shocks = zeros(S, nshocks, length(reg_ind))
        states[:, reg_ind], obs[:, reg_ind], pseudo[:, reg_ind], _ = forecast(system[reg_num], init_state, shocks)
    end

    # Use regime-switching forecast to iterate state-space system forward
    # during the forecast horizon
    if nperiods > regime_inds[end][end]
        fcast_shocks = zeros(S, nshocks, nperiods - regime_inds[end][end])
        fcast_inds = (regime_inds[end][end] + 1):nperiods

        states[:, fcast_inds], obs[:, fcast_inds], pseudo[:, fcast_inds], _ =
            forecast(m, system, states[:, regime_inds[end][end]], fcast_shocks;
                     cond_type = cond_type)
    end

    if start_index == 1 && end_index == nperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        return states[:, range], obs[:, range], pseudo[:, range]
    end
end

"""
```
trends(system::System{S}) where {S<:AbstractFloat}
trends(system::RegimeSwitchingSystem{S}) where {S<:AbstractFloat}
trends(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                start_date::Dates.Date = date_presample_start(m),
                end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                cond_type::Symbol) where {S<:AbstractFloat}
```

Compute trend (steady-state) states, observables, and pseudo-observables. The
trend is used for plotting shock decompositions. The first method
applies to non-regime switching systems, the second to regime-switching systems
that do not involve time variation in the `CCC` or `DD` vectors,
and the third to regime-switching systems with time variation in `CCC` or `DD`.
"""
function trends(system::System{S}) where {S<:AbstractFloat}

    state_trend  = system[:CCC]
    obs_trend    = system[:ZZ] * system[:CCC] + system[:DD]
    pseudo_trend = system[:ZZ_pseudo] * system[:CCC] + system[:DD_pseudo]

    return state_trend, obs_trend, pseudo_trend
end

function trends(system::RegimeSwitchingSystem{S}) where {S<:AbstractFloat}

    nreg  = n_regimes(system)
    input = 1:nreg

    state_trends  = Matrix{S}(undef, length(system[1, :CCC]), nreg)
    obs_trends    = Matrix{S}(undef, length(system[1, :DD]), nreg)
    pseudo_trends = Matrix{S}(undef, length(system[1, :DD_pseudo]), nreg)
    for i in input
        state_trends[:, i] .= 0.
        obs_trends[:, i]    = system[i, :DD]
        pseudo_trends[:, i] = system[i, :DD_pseudo]
    end

    return state_trends, obs_trends, pseudo_trends
end

# Handles case of time-varying CCC and DD
function trends(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                start_date::Dates.Date = date_presample_start(m),
                end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                cond_type::Symbol = :none) where {S<:AbstractFloat}

    # Dates: We compute the trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods    = subtract_quarters(date_forecast_end(m), date_mainsample_start(m)) + 1
    start_index = index_shockdec_start(m)
    end_index   = index_shockdec_end(m)

    # Now calculate the regime indices in history and forecast
    hist_regime_inds = regime_indices(m, start_date, end_date) # do not need to account for ZLB split b/c shocks don't matter for trends
    if hist_regime_inds[1][1] < 1
        hist_regime_inds[1] = 1:hist_regime_inds[1][end]
    end
    if hist_regime_inds[end][end] >= end_index  # if the end index is in the middle of a regime or is past the regime's end
        hist_regime_inds[end] = hist_regime_inds[end][1]:end_index
        fcast_regime_inds = nothing
    else
        fcast_regime_inds = get_fcast_regime_inds(m, forecast_horizons(m; cond_type = cond_type), cond_type,
                                                  start_index = hist_regime_inds[end][end])
        fcast_cutoff = findfirst([regind[end] >= end_index for regind in fcast_regime_inds])
        if isnothing(fcast_cutoff)
            error("The index_shockdec_end(m) occurs past the index of the final forecast date.")
        end
        fcast_regime_inds = fcast_regime_inds[1:fcast_cutoff]
        fcast_regime_inds[end] = fcast_regime_inds[end][1]:end_index
    end
    if length(hist_regime_inds[end]) == 0
        pop!(hist_regime_inds)
    end

    # Initialize output
    state_trends  = Matrix{S}(undef, length(system[1, :CCC]), end_index)
    obs_trends    = Matrix{S}(undef, length(system[1, :DD]), end_index)
    pseudo_trends = Matrix{S}(undef, length(system[1, :DD_pseudo]), end_index)

    # Check if any CCC are nonzero, otherwise we can do a faster version
    first_nonzero_CCC = findfirst([any(.!(system[i, :CCC] .≈ 0.)) for i in 1:n_regimes(system)])
    if isnothing(first_nonzero_CCC)
        state_trends .= 0.
        for (reg, inds) in enumerate(hist_regime_inds)
            for t in inds
                obs_trends[:, t]    = system[reg, :DD]
                pseudo_trends[:, t] = system[reg, :DD_pseudo]
            end
        end
        if !isnothing(fcast_regime_inds)
            regbuffer = cond_type == :none ? get_setting(m, :reg_forecast_start) - 1 :
                get_setting(m, :reg_post_conditional_end) - 1
            for (i, inds) in enumerate(fcast_regime_inds)
                reg = i + regbuffer
                for t in inds
                    obs_trends[:, t]    = system[reg, :DD]
                    pseudo_trends[:, t] = system[reg, :DD_pseudo]
                end
            end
        end
    else
        first_nonzero_CCC_hist, first_nonzero_CCC_fcast = if first_nonzero_CCC > length(hist_regime_inds)
            length(hist_regime_inds) + 1, first_nonzero_CCC - length(hist_regime_inds)
        else
            first_nonzero_CCC, 1
        end

        for (reg, inds) in enumerate(hist_regime_inds[1:(first_nonzero_CCC_hist - 1)])
            state_trends[:, inds]  .= 0.
            obs_trends[:, inds]    .= system[reg, :DD]
            pseudo_trends[:, inds] .= system[reg, :DD_pseudo]
        end

        for (i, inds) in enumerate(hist_regime_inds[first_nonzero_CCC_hist:end])
            reg = i + first_nonzero_CCC_hist - 1
            for t in inds
                state_trends[:, t]  = system[reg, :CCC] + system[reg, :TTT] * state_trends[:, t - 1]
                obs_trends[:, t]    = system[reg, :ZZ] * state_trends[:, t] + system[reg, :DD]
                pseudo_trends[:, t] = system[reg, :ZZ_pseudo] * state_trends[:, t] + system[reg, :DD_pseudo]
            end
        end

        if !isnothing(fcast_regime_inds)
            reg_buffer = cond_type == :none ? get_setting(m, :reg_forecast_start) - 1 :
                get_setting(m, :reg_post_conditional_end) - 1
            for (i, inds) in enumerate(fcast_regime_inds[1:(first_nonzero_CCC_fcast - 1)])
                reg = i + reg_buffer
                for t in inds
                    state_trends[:, inds]  .= 0.
                    obs_trends[:, inds]    .= system[reg, :DD]
                    pseudo_trends[:, inds] .= system[reg, :DD_pseudo]
                end
            end

            for (i, inds) in enumerate(fcast_regime_inds[first_nonzero_CCC_fcast:end])
                reg = i + first_nonzero_CCC_fcast - 1 + reg_buffer
                for t in inds
                    state_trends[:, t]  = system[reg, :CCC] + system[reg, :TTT] * state_trends[:, t - 1]
                    obs_trends[:, t]    = system[reg, :ZZ] * state_trends[:, t] + system[reg, :DD]
                    pseudo_trends[:, t] = system[reg, :ZZ_pseudo] * state_trends[:, t] + system[reg, :DD_pseudo]
                end
            end
        end
    end

    if start_index == 1
        return state_trends, obs_trends, pseudo_trends
    else
        return state_trends[:, start_index:end], obs_trends[:, start_index:end], pseudo_trends[:, start_index:end]
    end
end
