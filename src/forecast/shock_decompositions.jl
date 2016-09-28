"""
```
shock_decompositions{T<:AbstractFloat}(m::AbstractModel,
    syses::Vector{System{S}}, histshocks::Vector{Matrix{S}})
```

Computes shock decompositions for all draws, given a model object, system matrices, and
historical smoothed shocks.

### Inputs

- `m`: model object
- `syses::Vector{System}`: vector of length `ndraws`, whose elements are
  `System` objects specifying state-space system matrices for each draw
- `histshocks`::Vector{Matrix{S}}`: vector of length `ndraws`, whose elements
  are the `nshocks` x `hist_periods` matrices of historical smoothed shocks

### Outputs

-`states`: vector of length `ndraws`, whose elements are the `nstates` x
 `nperiods` x `nshocks` matrices of state shock decompositions
-`observables`: vector of length `ndraws`, whose elements are the `nobs` x
 `nperiods` x `nshocks` matrices of observable shock decompositions
-`pseudo_observables`: vector of length `ndraws`, whose elements are the
 `npseudo` x `nperiods` x `nshocks` matrices of pseudo-observable shock
 decompositions

where `nperiods = hist_periods + forecast_horizon`.
"""
function shock_decompositions{T<:AbstractFloat}(m::AbstractModel,
    syses::Vector{System{T}}, histshocks::Vector{Matrix{T}})

    ndraws = length(syses)

    # for now, we are ignoring pseudo-observables so these can be empty
    Z_pseudo = zeros(T, 12, n_states_augmented(m))
    D_pseudo = zeros(T, 12)

    # retrieve settings for forecast
    horizon = forecast_horizons(m)
    nshocks = n_shocks_exogenous(m)

    # Unpack everything for call to map/pmap
    TTTs     = map(s -> s[:TTT], syses)
    RRRs     = map(s -> s[:RRR], syses)
    ZZs      = map(s -> s[:ZZ],  syses)
    DDs      = map(s -> s[:DD],  syses)

    # Prepare copies of these objects due to lack of parallel broadcast functionality
    ZZps     = fill(Z_pseudo, ndraws)
    DDps     = fill(D_pseudo, ndraws)
    horizons = fill(horizon,  ndraws)

    # Determine periods for which to return shock decompositions
    starts = if !isnull(shockdec_startdate(m))
        start_index = subtract_quarters(get(shockdec_startdate(m)), date_prezlb_start(m)) + 1
        fill(start_index, ndraws)
    else
        fill(1, ndraws)
    end

    ends = if !isnull(shockdec_enddate(m))
        end_index = subtract_quarters(get(shockdec_enddate(m)), date_prezlb_start(m)) + 1
        fill(end_index, ndraws)
    else
        map(draw -> size(draw, 2) + horizon, histshocks)
    end

    # Go to work!
    if use_parallel_workers(m)
        mapfcn = pmap
    else
        mapfcn = map
    end

    shockdecs = mapfcn(DSGE.compute_shock_decompositions, TTTs, RRRs, ZZs, DDs,
                       ZZps, DDps, horizons, histshocks, starts, ends)

    # Unpack returned vector of tuples
    states             = [shockdec[1]::Array{T, 3} for shockdec in shockdecs]
    observables        = [shockdec[2]::Array{T, 3} for shockdec in shockdecs]
    pseudo_observables = [shockdec[3]::Array{T, 3} for shockdec in shockdecs]

    return states, observables, pseudo_observables
end


"""
```
compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S},
    Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
    forecast_horizons::Int, histshocks::Matrix{S}, start_index::Nullable{Int},
    end_index::Nullable{Int})
```

### Inputs

- `T`, `R`: transition equation matrices
- `Z`, `D`: observation equation matrices
- `Z_pseudo`, `D_pseudo`: matrices mapping states to pseudo-observables
- `forecast_horizons`: number of quarters ahead to forecast output
- `histshocks`: matrix of smoothed historical shocks (size `nshocks` x
  `hist_periods`)
- `start_index`: indicates the first index from which to return computed shock
  decompositions
- `end_index`: incidates the last index for which to return computed shock
  decompositions

### Outputs

`compute_shock_decompositions` returns a 3-tuple of shock decompositions:

-`:states`: `nstates` x `nperiods` x `nshocks`
-`:observables`: `nobservables` x `nperiods` x `nshocks`
-`:pseudo_observables`: `npseudo` x `nperiods` x `nshocks`

where `nperiods` is `end_index - start_index + 1`.
"""
function compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S},
    R::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, forecast_horizons::Int, histshocks::Matrix{S},
    start_index::Int, end_index::Int)

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobservables = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    histperiods  = size(histshocks, 2)
    allperiods   = histperiods + forecast_horizons

    if forecast_horizons <= 0 || start_index < 1 || end_index > allperiods
        throw(DomainError())
    end

    states             = zeros(S, nstates,      allperiods, nshocks)
    observables        = zeros(S, nobservables, allperiods, nshocks)
    pseudo_observables = zeros(S, npseudo,      allperiods, nshocks)

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

        # Apply observation and pseudo-observation equations
        observables[:, :, i]        = D        .+ Z        * states[:, :, i]
        pseudo_observables[:, :, i] = D_pseudo .+ Z_pseudo * states[:, :, i]
    end

    # Return shock decompositions in appropriate range
    if start_index == 1 && end_index == allperiods
        return states, observables, pseudo_observables
    else
        range = start_index:end_index
        return states[:, range, :], observables[:, range, :], pseudo_observables[:, range, :]
    end
end
