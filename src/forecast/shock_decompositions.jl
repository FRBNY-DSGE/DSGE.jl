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
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m));
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

    shock_decompositions(m, system, horizon, histshocks, start_index, end_index, regime_inds,
                         cond_type = cond_type)
end

function shock_decompositions(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S},
                              forecast_horizons::Int, histshocks::Matrix{S},
                              start_index::Int, end_index::Int,
                              regime_inds::Vector{UnitRange{Int}};
                              cond_type::Symbol = :none) where {S<:AbstractFloat}

    # forecast_horizons = 1
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

# Option 1: Passing in the model object
#=
        # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
        init_state = zeros(S, nstates) # update initial state

        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(m, system, init_state, shocks)
=#
# Option 2: Looping through each regime
        for (reg_num, reg_ind) in enumerate(regime_inds)
            if reg_ind[end] > histperiods
                break
            end
            init_state = (reg_num == 1) ? zeros(S, nstates) : states[:, reg_ind[1] - 1, i] # update initial state

            states[:, reg_ind, i], obs[:, reg_ind, i], pseudo[:, reg_ind, i], _ =
                forecast(system[reg_num], init_state, shocks[:, reg_ind])
        end

        if allperiods > histperiods
            fcast_inds = (histperiods + 1):allperiods
            states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
                forecast(m, system, states[:, histperiods, i], fcast_shocks, cond_type = cond_type)
        end

# Old Forecast setup
#=
        fcast_inds = (histperiods + 1):allperiods
        states[:, fcast_inds, i], obs[:, fcast_inds, i], pseudo[:, fcast_inds, i], _ =
            forecast(system[length(regime_inds)], states[:, regime_inds[end][end], i], fcast_shocks)
=#
        # zero out shocks
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
                              end_date::Dates.Date = prev_quarter(date_forecast_start(m));
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
    regime_inds[end] = regime_inds[end][1]:end_index # if the end index is in the middle of a regime or is past the regime's end
    if length(regime_inds[end]) == 0
        pop!(regime_inds)
    end

    return deterministic_trends(m, system, z0, nperiods, start_index, end_index, regime_inds, cond_type = cond_type)
end


function deterministic_trends(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S}, z0::Vector{S}, nperiods::Int,
                              start_index::Int, end_index::Int, regime_inds::Vector{UnitRange{Int}};
                              cond_type::Symbol = :none) where {S<:AbstractFloat}

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
#=
# Option 1: Passing in model object

    shocks = zeros(S, nshocks, nperiods)
    states, obs, pseudo, _ = forecast(m, system, z0, shocks)
=#

# Option 2: Forecasting for each regime via loop
    # Use forecast to iterate state-space system forward without shocks or the ZLB procedure
    for (reg_num, reg_ind) in enumerate(regime_inds)
        if reg_num <= get_setting(m, :n_hist_regimes)
            init_state = (reg_num == 1) ? z0 : states[:, reg_ind[1] - 1] # update initial state

            # Construct matrix of 0 shocks for this regime
            shocks = zeros(S, nshocks, length(reg_ind))
            states[:, reg_ind], obs[:, reg_ind], pseudo[:, reg_ind], _ = forecast(system[reg_num], init_state, shocks)
        end
    end

    # fcast_inds = (get_setting(m, :n_hist_regimes) + 1):get_setting(m, :n_regimes)
    if nperiods > regime_inds[end][end]
        fcast_shocks = zeros(S, nshocks, nperiods - regime_inds[end][end])
        fcast_inds = (regime_inds[end][end] + 1):nperiods
        states[:, fcast_inds], obs[:, fcast_inds], pseudo[:, fcast_inds], _ =
            forecast(m, system, states[:, regime_inds[end][end]], fcast_shocks;
                     cond_type = cond_type)
    end
#=
    if regime_inds[end][end] < nperiods
        shocks   = zeros(S, nshocks, nperiods - regime_inds[end][end])
        last_ind = (1 + regime_inds[end][end]):nperiods
        states[:, last_ind], obs[:, last_ind], pseudo[:, last_ind], _ = forecast(system[length(regime_inds)], states[:, regime_inds[end][end]], shocks)
    end
=#
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

function trends(system::RegimeSwitchingSystem{S}, m, start_date, end_date) where {S<:AbstractFloat}

    nreg  = n_regimes(system)
    input = 1:nreg

    dated = start_date
    nper_vec = zeros(Int64, n_regimes(system))
    for i in 1:n_regimes(system)
        dated = get_setting(m, :regime_dates)[i]
        if i != n_regimes(system)
            dated2 = get_setting(m, :regime_dates)[i+1]
            nper_vec[i] = DSGE.subtract_quarters(dated2, dated)
        else
            nper_vec[i] = max(DSGE.subtract_quarters(Date(2025,12,31), dated) + 1, 1)
        end
    end

    state_trends  = Matrix{S}(undef, length(system[1, :CCC]), sum(nper_vec))
    obs_trends    = Matrix{S}(undef, length(system[1, :DD]), sum(nper_vec))
    pseudo_trends = Matrix{S}(undef, length(system[1, :DD_pseudo]), sum(nper_vec))

    @show nper_vec
    for j in 1:length(nper_vec)
        for k in 1:nper_vec[j]
            if j == 1 && k == 1
                state_trends[:, 1]  = system[j, :CCC]
                obs_trends[:, 1]    = system[j, :ZZ] * system[j, :CCC] + system[j, :DD]
                pseudo_trends[:, 1] = system[j, :ZZ_pseudo] * system[j, :CCC] + system[j, :DD_pseudo]
            else
                counter = j == 1 ? k : sum(nper_vec[1:(j-1)]) + k
                @show j, k
                state_trends[:, counter]  = system[j, :CCC] + system[j, :TTT] * state_trends[:,counter-1]
                obs_trends[:, counter]    = system[j, :ZZ] * state_trends[:,counter] + system[j, :DD]
                pseudo_trends[:, counter] = system[j, :ZZ_pseudo] * state_trends[:,counter] + system[j, :DD_pseudo]
            end
        end
    end

    return state_trends, obs_trends, pseudo_trends
end
