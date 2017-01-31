"""
```
compute_shock_decompositions(system, forecast_horizons, histshocks, start_index,
    end_index)

compute_shock_decompositions(T, R, Z, D, Z_pseudo, D_pseudo, z0, shocks,
    forecast_horizons, histshocks, start_index, end_index)
```

### Inputs

- `system::System{S}`: state-space system matrices. Alternatively, provide
  transition equation matrices `T`, `R`; measurement equation matrices `Z`, `D`;
  and (possibly empty) pseudo-measurement equation matrices `Z_pseudo` and
  `D_pseudo`.
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
  pseudo-observable shock decompositions. If the provided `Z_pseudo` and
  `D_pseudo` matrices are empty, then `pseudo` will be empty.

where `nperiods = `end_index - start_index + 1`.
"""
function compute_shock_decompositions{S<:AbstractFloat}(system::System{S},
    forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int)

    # Unpack system
    T, R = system[:TTT], system[:RRR]
    Z, D = system[:ZZ], system[:DD]

    Z_pseudo, D_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    compute_shock_decompositions(T, R, Z, D, Z_pseudo, D_pseudo,
        forecast_horizons, histshocks, start_index, end_index)
end

function compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S},
    R::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int)

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobs         = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    histperiods  = size(histshocks, 2)
    allperiods   = histperiods + forecast_horizons

    forecast_pseudo = !isempty(Z_pseudo) && !isempty(D_pseudo)

    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError())
    end

    states = zeros(S, nstates,      allperiods, nshocks)
    obs    = zeros(S, nobs, allperiods, nshocks)
    pseudo = if forecast_pseudo
        zeros(S, npseudo, allperiods, nshocks)
    else
        zeros(S, 0, 0, 0)
    end

    # Define our iteration function
    iterate(z_t1, ϵ_t) = T*z_t1 + R*ϵ_t
    z0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, allperiods)
        shocks[i, 1:histperiods] = histshocks[i, :]

        # Iterate state space forward
        states[:, 1, i] = iterate(z0, shocks[:, 1])
        for t in 2:allperiods
            states[:, t, i] = iterate(states[:, t-1, i], shocks[:, t])
        end

        # Apply measurement and pseudo-measurement equations
        obs[:, :, i] = D .+ Z * states[:, :, i]
        if forecast_pseudo
            pseudo[:, :, i] = D_pseudo .+ Z_pseudo * states[:, :, i]
        end

    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, obs, pseudo
    else
        range = start_index:end_index
        if forecast_pseudo
            return states[:, range, :], obs[:, range, :], pseudo[:, range, :]
        else
            return states[:, range, :], obs[:, range, :], pseudo
        end
    end
end

"""
```
compute_deterministic_trend(m, system, z0)

compute_deterministic_trend(system, z0, nperiods)
```

Compute deterministic trend values of states, observables, and
pseudoobservables, given a model object and system matrices. The deterministic
trend for a single draw is simply the series that would be obtained by iterating
the state-space system forward, beginning from a state vector `z0` in the last
presample period.

### Inputs

- `m::AbstractModel`: model object
- `system::System{S}`: state-space system matrices
- `z0`::Vector{S}: initial state vector

where `S<:AbstractFloat`.

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `nperiods` of state
  steady-state values
- `obs::Matrix{S}`: matrix of size `nobs` x `nperiods` of observable
  steady-state values
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `nperiods` of
  pseudo-observable steady-state values. If `!forecast_pseudoobservables(m)`,
  `pseudo` will be empty.

where `nperiods` is the number of quarters between
`date_mainsample_start(m)` and `date_shockdec_end(m)`, inclusive.
"""
function deterministic_trends{S<:AbstractFloat}(m::AbstractModel{S},
    system::System{S}, z0::Vector{S})

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods = subtract_quarters(date_shockdec_end(m), date_mainsample_start(m)) + 1

    compute_deterministic_trend(system, z0, nperiods)
end

function compute_deterministic_trend{S<:AbstractFloat}(system::System{S}, z0::Vector{S}, nperiods::Int)

    # construct matrix of 0 shocks for entire history and forecast horizon
    nshocks  = size(system[:RRR], 2)
    shocks   = zeros(S, nshocks, nperiods)

    # use compute_forecast to iterate state-space system forward without shocks or the ZLB procedure
    states, obs, pseudo, _ = compute_forecast(system, z0, shocks)

    states, obs, pseudo
end

"""
```
compute_trend{S<:AbstractFloat}(system::System{S})
```

Compute trend (steady-state) states, observables, and pseudo-observables. The
trend is used for plotting shock decompositions.
"""
function compute_trend{S<:AbstractFloat}(system::System{S})

    # Unpack system
    C, D = system[:CCC], system[:DD]

    D_pseudo = if !isnull(system.pseudo_measurement)
        system[:DD_pseudo]
    else
        Vector{S}()
    end

    C, D, D_pseudo
end