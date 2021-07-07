"""
```
shock_decompositions(m, system, histshocks)

shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)

shock_decompositions(m, system, histshocks, start_date, end_date)

shock_decompositions(m, system, forecast_horizons, histshocks, start_index,
    end_index, regime_inds, cond_type)

shock_decompositions(m, system, old_system,
                              histshocks, old_histshocks, start_date, end_date,
                              cond_type; full_shock_decomp)

shock_decompositions(m, system, old_system,
                              forecast_horizons, histshocks, old_histshocks,
                              start_index, end_index,
                              regime_inds, cond_type; full_shock_decomp)
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
- `old_system::RegimeSwitchingSystem{S}`: state-space system matrices for the old model
- `old_histshocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of
  historical smoothed shocks using the old system
- `full_shock_decomp::Bool`: If true for old_system case, return difference in shockdecs
    between all shocks for the new and old systems. Else, return old system's shocks
    with 1 cumulative AIT shock added on.

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
## effect of using new system compared to old.
function shock_decompositions(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, old_system::RegimeSwitchingSystem{S},
                              histshocks::Matrix{S}, old_histshocks::Matrix{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none; full_shock_decomp = true) where {S<:AbstractFloat}
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

    shock_decompositions(m, system, old_system, horizon, histshocks, old_histshocks, start_index, end_index,
                         regime_inds, cond_type, full_shock_decomp = full_shock_decomp)
end

function shock_decompositions(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                              old_system::RegimeSwitchingSystem{S},
                              forecast_horizons::Int,
                              histshocks::Matrix{S}, old_histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol; full_shock_decomp = true) where {S<:AbstractFloat}

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
    old_shocks = zeros(S, nshocks, histperiods)

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
    old_system = zero_system_constants(old_system)
    z0     = zeros(S, nstates)

    fcast_shocks = zeros(S, nshocks, forecast_horizons) # these are always zero
    for i = 1:nshocks
        # Isolate single shock
        shocks[i, :] = histshocks[i, :]
        old_shocks[i, :] = old_histshocks[i, :]

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
            #=if reg_ind[end] < old_system
                old_system2 = system[reg_num]
                counts = reg_num
            else
                old_system2 = system[1]#system[counts]
            end=#
            old_system2 = old_system
            init_state_old = (reg_num == 1) ? zeros(S, nstates) : old_states[:, reg_ind[1] - 1, i]
            old_states[:, reg_ind, i], old_obs[:, reg_ind, i], old_pseudo[:, reg_ind, i], _ =
                forecast(old_system2[reg_num], init_state_old, old_shocks[:,reg_ind])
        end

        # Run forecast using regime-switching forecast
        if allperiods > histperiods
            fcast_inds = (histperiods + 1):allperiods
            states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
                forecast(m, system, states[:, histperiods, i], fcast_shocks, cond_type = cond_type)

            # Old System Values
            old_system2 = old_system#system[1]#system[counts]
            old_states[:, fcast_inds, i], old_obs[:, fcast_inds, i], old_pseudo[:, fcast_inds, i], _ =
                forecast(m, old_system2, old_states[:, histperiods, i], fcast_shocks, cond_type = cond_type)
                #forecast(old_system2, old_states[:, histperiods, i], fcast_shocks)
        end

        # zero out shocks b/c want effects of single shocks
        shocks[i, :] .= 0.
        old_shocks[i, :] .= 0.
    end

    # Add new system as shock to old system
    #=old_states = cat(old_states, zeros(size(old_states[:,:,1])), dims = 3)
    old_obs = cat(old_obs, zeros(size(old_obs[:,:,1])), dims = 3)
    old_pseudo = cat(old_pseudo, zeros(size(old_pseudo[:,:,1])), dims = 3)=#

    for i in 1:size(states,1)
        for j in 1:size(states,2)
            old_states[i,j,end] = sum([states[i,j,k] for k in 1:(size(states,3)-1)]) - sum([old_states[i,j,k] for k in 1:(size(states,3)-1)])
            full_shock_decomp && (old_states[i,j,1:end-1] .= states[i,j,1:end-1] .- old_states[i,j,1:end-1])
        end
    end

    for i in 1:size(obs,1)
        for j in 1:size(obs,2)
            old_obs[i,j,end] = sum([obs[i,j,k] for k in 1:(size(obs,3)-1)]) - sum([old_obs[i,j,k] for k in 1:(size(obs,3)-1)])
            full_shock_decomp && (old_obs[i,j,1:end-1] .= obs[i,j,1:end-1] .- old_obs[i,j,1:end-1])
        end
    end

    for i in 1:size(pseudo,1)
        for j in 1:size(pseudo,2)
            old_pseudo[i,j,end] = sum([pseudo[i,j,k] for k in 1:(size(pseudo,3)-1)]) - sum([old_pseudo[i,j,k] for k in 1:(size(pseudo,3)-1)])
            full_shock_decomp && (old_pseudo[i,j,1:end-1] .= pseudo[i,j,1:end-1] .- old_pseudo[i,j,1:end-1])
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

deterministic_trends(m, system, old_system, z0, start_date, end_date)

deterministic_trends(m, system, old_system, z0, nperiods, start_index, end_index,
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
- `old_system::`RegimeSwitchingSystem`: state-space system matrices for old model

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

# Get Deterministic Trends with an old system
function deterministic_trends(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, old_system::RegimeSwitchingSystem{S},
                              z0::Vector{S},
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

    return deterministic_trends(m, system, old_system, z0, nperiods, start_index, end_index, regime_inds, cond_type)
end


function deterministic_trends(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S}, old_system::RegimeSwitchingSystem{S},
                              z0::Vector{S}, nperiods::Int,
                              start_index::Int, end_index::Int, regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol) where {S<:AbstractFloat}

    # Set constant system matrices to 0
    system = zero_system_constants(system)
    old_system = zero_system_constants(old_system)

    # Initialize storage matrices
    nshocks = size(system[1, :RRR], 2) # this is used in the forecasting loops to create shock matrices
    nstates = size(system[1, :TTT], 1)
    nobs    = size(system[1, :ZZ], 1)
    npseudo = size(system[1, :ZZ_pseudo], 1)
    states = Matrix{S}(undef, nstates, nperiods)
    obs    = Matrix{S}(undef, nobs,    nperiods)
    pseudo = Matrix{S}(undef, npseudo, nperiods)

    # Old System
    states_old = Matrix{S}(undef, nstates, nperiods)
    obs_old    = Matrix{S}(undef, nobs,    nperiods)
    pseudo_old = Matrix{S}(undef, npseudo, nperiods)

    # Use forecast to iterate state-space system forward
    # without shocks or the ZLB procedure during history
    for (reg_num, reg_ind) in enumerate(regime_inds)
        init_state = (reg_num == 1) ? z0 : states[:, reg_ind[1] - 1] # update initial state

        # Construct matrix of 0 shocks for this regime
        shocks = zeros(S, nshocks, length(reg_ind))
        states[:, reg_ind], obs[:, reg_ind], pseudo[:, reg_ind], _ = forecast(system[reg_num], init_state, shocks)

        # Do for old system too
        old_system2 = old_system[reg_num]#reg_ind[end] >= old_system ? system[1] : system[reg_num]
        init_state_old = (reg_num == 1) ? z0 : states_old[:, reg_ind[1] - 1] # update initial state
        states_old[:, reg_ind], obs_old[:, reg_ind], pseudo_old[:, reg_ind], _ = forecast(old_system2, init_state, shocks)
    end

    # Use regime-switching forecast to iterate state-space system forward
    # during the forecast horizon
    if nperiods > regime_inds[end][end]
        fcast_shocks = zeros(S, nshocks, nperiods - regime_inds[end][end])
        fcast_inds = (regime_inds[end][end] + 1):nperiods

        states[:, fcast_inds], obs[:, fcast_inds], pseudo[:, fcast_inds], _ =
            forecast(m, system, states[:, regime_inds[end][end]], fcast_shocks;
                     cond_type = cond_type)

        # Old System
        ##TODO: (Same in shock_decompositions) - currently assuming that we use
        ## Taylor in all forecast regimes, but generalize based on regime_inds[end][end]
        ## and old_system values
        states_old[:, fcast_inds], obs_old[:, fcast_inds], pseudo_old[:, fcast_inds], _ =
            forecast(m, old_system, states_old[:, regime_inds[end][end]], fcast_shocks, cond_type = cond_type)
        #forecast(old_system, states_old[:, regime_inds[end][end]], fcast_shocks)
    end

    # Old System Update
    states = states .- states_old
    obs = obs .- obs_old
    pseudo = pseudo .- pseudo_old

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

# Sequence of Shock Decompositions (Effect of time t shock on obs/pseudo/states at each time for all t)
"""
```
shock_decompositions_sequence(m, system, histshocks, n_back, back_shocks)

shock_decompositions_sequence(system, forecast_horizons, histshocks, start_index,
    end_index, n_back, back_shocks)

shock_decompositions_sequence(m, system, histshocks, start_date, end_date, n_back, back_shocks)

shock_decompositions_sequence(m, system, forecast_horizons, histshocks, start_index,
    end_index, regime_inds, cond_type, n_back, back_shocks)

shock_decompositions_sequence(m, system, old_system,
                              histshocks, old_histshocks, start_date, end_date,
                              cond_type; full_shock_decomp, n_back, back_shocks)

shock_decompositions_sequence(m, system, old_system,
                              forecast_horizons, histshocks, old_histshocks,
                              start_index, end_index,
                              regime_inds, cond_type; full_shock_decomp, n_back, back_shocks)
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
- `old_system::RegimeSwitchingSystem{S}`: state-space system matrices for the old model
- `old_histshocks::Matrix{S}`: matrix of size `nshocks` x `hist_periods` of
  historical smoothed shocks using the old system
- `full_shock_decomp::Bool`: If true for old_system case, return difference in shockdecs
    between all shocks for the new and old systems. Else, return old system's shocks
    with 1 cumulative AIT shock added on.

where `S<:AbstractFloat`.

### Outputs

- `states::Array{S, 4}`: matrix of size `nhistperiods` x `nstates` x `nperiods` x `nshocks` of
  state shock decompositions
- `obs::Array{S, 4}`: matrix of size `nhistperiods` x `nobs` x `nperiods` x `nshocks` of
  observable shock decompositions
- `pseudo::Array{S, 4}`: matrix of size `nhistperiods` x `npseudo` x `nperiods` x `nshocks` of
  pseudo-observable shock decompositions

where `nperiods = `end_index - start_index + 1`
and `nhistperiods` = `size(histshocks,2)`
"""
function shock_decompositions_sequence(m::AbstractDSGEModel,
    system::System{S}, histshocks::Matrix{S};
    n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}

    horizon   = forecast_horizons(m)
    start_ind = index_shockdec_start(m)
    end_ind   = index_shockdec_end(m)

    shock_decompositions_sequence(system, horizon, histshocks, start_ind, end_ind, n_back = n_back, back_shocks = back_shocks)
end

function shock_decompositions_sequence(system::System{S},
    forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int;
    n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}

    # Setup
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 2)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)
    histperiods  = size(histshocks, 2)
    allperiods   = histperiods + forecast_horizons

    dim3_len = min(end_index-start_index+1, allperiods)
    states = zeros(S, histperiods, nstates, dim3_len, nshocks)
    obs    = zeros(S, histperiods, nobs, dim3_len, nshocks)
    pseudo = zeros(S, histperiods, npseudo, dim3_len, nshocks)

    for t in 1:histperiods
        histshocks_t = copy(histshocks)
        histshocks_t[:,vcat(1:t-1,t+1:histperiods)] .= 0.0

        states[t,:,:,:], obs[t,:,:,:], pseudo[t,:,:,:] = shock_decompositions(system, forecast_horizons,
                                                                                  histshocks_t, start_index, end_index)
    end

    # Get the matrices back to 3 dimensions by having bars for contemporaneous shocks and shocks each of n_back time periods back.
    ## n_back is given by user and shocks for which we don't just have contemporaneous shocks is given by user.
    states3 = zeros(size(states,2), size(states,3), size(states,4) + (n_back+1)*length(back_shocks))
    obs3 = zeros(size(obs,2), size(obs,3), size(obs,4) + (n_back+1)*length(back_shocks))
    pseudo3 = zeros(size(pseudo,2), size(pseudo,3), size(pseudo,4) + (n_back+1)*length(back_shocks))

    ## Get indices for vector of shocks
    shock_inds = [m.exogenous_shocks[shock_ind] for shock_ind in back_shocks]
    shock_len = length(back_shocks)

    ## Actually do the conversion
    for i in 1:size(states,3)
        st_ind = min(start_index+i-1, size(states,1))
        states3[:,i,1:size(states,4)] = states[st_ind,:,i,:]
        obs3[:,i,1:size(states,4)] = obs[st_ind,:,i,:]
        pseudo3[:,i,1:size(states,4)] = pseudo[st_ind,:,i,:]

        for k in 1:n_back
            st_ind2 = min(start_index+i-1-k, size(states,1)-k)
            for j in 1:length(back_shocks)
                states3[:,i,size(states,4)+(k-1)*shock_len+j] = states[st_ind2,:,i,shock_inds[j]]
                obs3[:,i,size(obs,4)+(k-1)*shock_len+j] = obs[st_ind2,:,i,shock_inds[j]]
                pseudo3[:,i,size(pseudo,4)+(k-1)*shock_len+j] = pseudo[st_ind2,:,i,shock_inds[j]]

                # Add this "lagged shock" to the list of exogenous shocks
                if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_, k))
                    m.exogenous_shocks[Symbol(back_shocks[j], :_, k)] = length(m.exogenous_shocks)+1
                end
            end
        end
        ## Add in a bar that sums up shocks before n_back periods back.
        st_ind3 = min(start_index+i-2-n_back, size(states,1)-n_back-1)
        for j in 1:length(back_shocks)
            states3[:,i,size(states,4)+n_back*shock_len+j] = sum(states[1:st_ind3,:,i,shock_inds[j]],dims=1)
            obs3[:,i,size(obs,4)+n_back*shock_len+j] = sum(obs[1:st_ind3,:,i,shock_inds[j]],dims=1)
            pseudo3[:,i,size(pseudo,4)+n_back*shock_len+j] = sum(pseudo[1:st_ind3,:,i,shock_inds[j]],dims=1)

            # Add this "lagged shock" to the list of exogenous shocks
            if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_History))
                m.exogenous_shocks[Symbol(back_shocks[j], :_History)] = length(m.exogenous_shocks)+1
            end
        end
    end

    return states3, obs3, pseudo3
end


function shock_decompositions_sequence(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, histshocks::Matrix{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none;
                              n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}
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

    shock_decompositions_sequence(m, system, horizon, histshocks, start_index, end_index,
                                  regime_inds, cond_type, n_back = n_back, back_shocks = back_shocks)
end

function shock_decompositions_sequence(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                              forecast_horizons::Int, histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}},
                              cond_type::Symbol;
                              n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}

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

    dim3_len = min(end_index-start_index+1, allperiods)
    states = zeros(S, histperiods, nstates, dim3_len, nshocks)
    obs    = zeros(S, histperiods, nobs, dim3_len, nshocks)
    pseudo = zeros(S, histperiods, npseudo, dim3_len, nshocks)

    for t in 1:histperiods
        histshocks_t = copy(histshocks)
        histshocks_t[:,vcat(1:t-1,t+1:histperiods)] .= 0.0

        states[t,:,:,:], obs[t,:,:,:], pseudo[t,:,:,:] = shock_decompositions(m, system, forecast_horizons,
                                                                              histshocks_t, start_index, end_index,
                                                                              regime_inds, cond_type)
    end

    # Get the matrices back to 3 dimensions by having bars for contemporaneous shocks and shocks each of n_back time periods back.
    ## n_back is given by user and shocks for which we don't just have contemporaneous shocks is given by user.
    states3 = zeros(size(states,2), size(states,3), size(states,4) + (n_back+1)*length(back_shocks))
    obs3 = zeros(size(obs,2), size(obs,3), size(obs,4) + (n_back+1)*length(back_shocks))
    pseudo3 = zeros(size(pseudo,2), size(pseudo,3), size(pseudo,4) + (n_back+1)*length(back_shocks))

    ## Get indices for vector of shocks
    shock_inds = [m.exogenous_shocks[shock_ind] for shock_ind in back_shocks]
    shock_len = length(back_shocks)

    ## Actually do the conversion
    for i in 1:size(states,3)
        st_ind = min(start_index+i-1, size(states,1))
        states3[:,i,1:size(states,4)] = states[st_ind,:,i,:]
        obs3[:,i,1:size(states,4)] = obs[st_ind,:,i,:]
        pseudo3[:,i,1:size(states,4)] = pseudo[st_ind,:,i,:]

        for k in 1:n_back
            st_ind2 = min(start_index+i-1-k, size(states,1)-k)
            for j in 1:length(back_shocks)
                states3[:,i,size(states,4)+(k-1)*shock_len+j] = states[st_ind2,:,i,shock_inds[j]]
                obs3[:,i,size(obs,4)+(k-1)*shock_len+j] = obs[st_ind2,:,i,shock_inds[j]]
                pseudo3[:,i,size(pseudo,4)+(k-1)*shock_len+j] = pseudo[st_ind2,:,i,shock_inds[j]]

                # Add this "lagged shock" to the list of exogenous shocks
                if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_, k))
                    m.exogenous_shocks[Symbol(back_shocks[j], :_, k)] = length(m.exogenous_shocks)+1
                end
            end
        end
        ## Add in a bar that sums up shocks before n_back periods back.
        st_ind3 = min(start_index+i-2-n_back, size(states,1)-n_back-1)
        for j in 1:length(back_shocks)
            states3[:,i,size(states,4)+n_back*shock_len+j] = sum(states[1:st_ind3,:,i,shock_inds[j]],dims=1)
            obs3[:,i,size(obs,4)+n_back*shock_len+j] = sum(obs[1:st_ind3,:,i,shock_inds[j]],dims=1)
            pseudo3[:,i,size(pseudo,4)+n_back*shock_len+j] = sum(pseudo[1:st_ind3,:,i,shock_inds[j]],dims=1)

            # Add this "lagged shock" to the list of exogenous shocks
            if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_History))
                m.exogenous_shocks[Symbol(back_shocks[j], :_History)] = length(m.exogenous_shocks)+1
            end
        end
    end

    return states3, obs3, pseudo3
end

# Get Shock Decompositions with an old system and a bar that gives
## effect of using new system compared to old.
function shock_decompositions_sequence(m::AbstractDSGEModel{S},
                              system::RegimeSwitchingSystem{S}, old_system::RegimeSwitchingSystem{S},
                              histshocks::Matrix{S}, old_histshocks::Matrix{S},
                              start_date::Dates.Date = date_presample_start(m),
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m)),
                              cond_type::Symbol = :none; full_shock_decomp = true,
                              n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}
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

    shock_decompositions_sequence(m, system, old_system, horizon, histshocks, old_histshocks, start_index, end_index,
                                  regime_inds, cond_type, full_shock_decomp = full_shock_decomp,
                                  n_back = n_back, back_shocks = back_shocks)
end

function shock_decompositions_sequence(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                                       old_system::RegimeSwitchingSystem{S}, forecast_horizons::Int,
                                       histshocks::Matrix{S}, old_histshocks::Matrix{S},
                                       start_index::Int, end_index::Int,
                                       regime_inds::Vector{UnitRange{Int}},
                                       cond_type::Symbol; full_shock_decomp = true,
                                       n_back::Int64 = 0, back_shocks::Vector{Symbol} = Symbol[]) where {S<:AbstractFloat}

    # Setup
    nshocks     = size(system[1, :RRR], 2)
    nstates     = size(system[1, :TTT], 2)
    nobs        = size(system[1, :ZZ], 1)
    npseudo     = size(system[1, :ZZ_pseudo], 1)
    histperiods = size(histshocks, 2)
    allperiods  = histperiods + forecast_horizons

    dim3_len = min(end_index-start_index+1, allperiods)
    states = zeros(S, histperiods, nstates, dim3_len, nshocks)
    obs    = zeros(S, histperiods, nobs, dim3_len, nshocks)
    pseudo = zeros(S, histperiods, npseudo, dim3_len, nshocks)

    for t in 1:histperiods
        histshocks_t = copy(histshocks)
        histshocks_t[:,vcat(1:t-1,t+1:histperiods)] .= 0.0

        old_histshocks_t = copy(old_histshocks)
        old_histshocks_t[:,vcat(1:t-1,t+1:histperiods)] .= 0.0

        states[t,:,:,:], obs[t,:,:,:], pseudo[t,:,:,:] = shock_decompositions(m, system, old_system,
                                                                              forecast_horizons,
                                                                              histshocks_t, old_histshocks_t,
                                                                              start_index, end_index,
                                                                              regime_inds, cond_type,
                                                                              full_shock_decomp = full_shock_decomp)
    end

    # Get the matrices back to 3 dimensions by having bars for contemporaneous shocks and shocks each of n_back time periods back.
    ## n_back is given by user and shocks for which we don't just have contemporaneous shocks is given by user.
    states3 = zeros(size(states,2), size(states,3), size(states,4) + (n_back+1)*length(back_shocks))
    obs3 = zeros(size(obs,2), size(obs,3), size(obs,4) + (n_back+1)*length(back_shocks))
    pseudo3 = zeros(size(pseudo,2), size(pseudo,3), size(pseudo,4) + (n_back+1)*length(back_shocks))

    ## Get indices for vector of shocks
    shock_inds = [m.exogenous_shocks[shock_ind] for shock_ind in back_shocks]
    shock_len = length(back_shocks)

    ## Actually do the conversion
    for i in 1:size(states,3)
        st_ind = min(start_index+i-1, size(states,1))
        states3[:,i,1:size(states,4)] = states[st_ind,:,i,:]
        obs3[:,i,1:size(states,4)] = obs[st_ind,:,i,:]
        pseudo3[:,i,1:size(states,4)] = pseudo[st_ind,:,i,:]

        for k in 1:n_back
            st_ind2 = min(start_index+i-1-k, size(states,1)-k)
            for j in 1:length(back_shocks)
                states3[:,i,size(states,4)+(k-1)*shock_len+j] = states[st_ind2,:,i,shock_inds[j]]
                obs3[:,i,size(obs,4)+(k-1)*shock_len+j] = obs[st_ind2,:,i,shock_inds[j]]
                pseudo3[:,i,size(pseudo,4)+(k-1)*shock_len+j] = pseudo[st_ind2,:,i,shock_inds[j]]

                # Add this "lagged shock" to the list of exogenous shocks
                if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_, k))
                    m.exogenous_shocks[Symbol(back_shocks[j], :_, k)] = length(m.exogenous_shocks)+1
                end
            end
        end
        ## Add in a bar that sums up shocks before n_back periods back.
        st_ind3 = min(start_index+i-2-n_back, size(states,1)-n_back-1)
        for j in 1:length(back_shocks)
            states3[:,i,size(states,4)+n_back*shock_len+j] = sum(states[1:st_ind3,:,i,shock_inds[j]],dims=1)
            obs3[:,i,size(obs,4)+n_back*shock_len+j] = sum(obs[1:st_ind3,:,i,shock_inds[j]],dims=1)
            pseudo3[:,i,size(pseudo,4)+n_back*shock_len+j] = sum(pseudo[1:st_ind3,:,i,shock_inds[j]],dims=1)

            # Add this "lagged shock" to the list of exogenous shocks
            if !haskey(m.exogenous_shocks, Symbol(back_shocks[j], :_History))
                m.exogenous_shocks[Symbol(back_shocks[j], :_History)] = length(m.exogenous_shocks)+1
            end
        end
    end

    return states3, obs3, pseudo3
end
