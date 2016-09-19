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

- `states`: 4-dimensional array of size `nstates` x `nperiods` x `nshocks` x
 `ndraws` consisting of state shock decompositions for each draw
- `observables`: 4-dimensional array of size `nobs` x `nperiods` x `nshocks` x
 `ndraws` consisting of observable shock decompositions for each draw
- `pseudo_observables`: 4-dimensional array of size `npseudo` x `nperiods` x
 `nshocks` x `ndraws` consisting of pseudo-observable shock decompositions for
 each draw

where `nperiods = hist_periods + forecast_horizon`.
"""
function shock_decompositions{T<:AbstractFloat}(m::AbstractModel,
    syses::Vector{System{T}}, histshocks::Vector{Matrix{T}})

    ndraws = length(syses)

    # for now, we are ignoring pseudo-observables so these can be empty
    Z_pseudo = zeros(T, 12, n_states_augmented(m))
    D_pseudo = zeros(T, 12)

    # retrieve settings for forecast
    horizon  = forecast_horizons(m)
    nshocks  = n_shocks_exogenous(m)

    # Unpack everything for call to map/pmap
    TTTs   = [s[:TTT] for s in syses]
    RRRs   = [s[:RRR] for s in syses]
    ZZs    = [s[:ZZ]  for s in syses]
    DDs    = [s[:DD]  for s in syses]

    # Prepare copies of these objects due to lack of parallel broadcast functionality
    ZZps     = fill(Z_pseudo, ndraws)
    DDps     = fill(D_pseudo, ndraws)
    horizons = fill(horizon,  ndraws)

    if use_parallel_workers(m)
        mapfcn = pmap
    else
        mapfcn = map
    end

    # Go to work!
    shockdecs = mapfcn(DSGE.compute_shock_decompositions, TTTs, RRRs, ZZs, DDs,
        ZZps, DDps, horizons, histshocks)

    # Unpack returned vector of tuples
    states             = [shockdec[1]::Array{T, 3} for shockdec in shockdecs]
    observables        = [shockdec[2]::Array{T, 3} for shockdec in shockdecs]
    pseudo_observables = [shockdec[3]::Array{T, 3} for shockdec in shockdecs]

    # Splat vectors of 3-D arrays into 4-D arrays
    states             = cat(4, states...)
    observables        = cat(4, observables...)
    pseudo_observables = cat(4, pseudo_observables...)

    return states, observables, pseudo_observables
end


"""
```
compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S},
    C::Vector{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, forecast_horizons::Int, histshocks::Matrix{S})
```

### Inputs

- `T`, `R`: transition equation matrices
- `Z`, `D`: observation equation matrices
- `Z_pseudo`, `D_pseudo`: matrices mapping states to pseudo-observables
- `forecast_horizons`: number of quarters ahead to forecast output
- `histshocks`: matrix of smoothed historical shocks (size `nshocks` x
  `hist_periods`)

### Outputs

`compute_shock_decompositions` returns a 3-tuple of shock decompositions over
the both the history and the forecast period, whose elements have sizes:

- `:states`: `nstates` x `nperiods` x `nshocks`
- `:observables`: `nobs` x `nperiods` x `nshocks`
- `:pseudo_observables`: `npseudo` x `nperiods` x `nshocks`

where `nperiods = hist_periods + forecast_horizon`.
"""
function compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S},
    R::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, forecast_horizons::Int, histshocks::Matrix{S})

    if forecast_horizons <= 0
        throw(DomainError())
    end

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobservables = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    histperiods  = size(histshocks, 2)
    nperiods     = histperiods + forecast_horizons

    states             = zeros(S, nstates,      nperiods, nshocks)
    observables        = zeros(S, nobservables, nperiods, nshocks)
    pseudo_observables = zeros(S, npseudo,      nperiods, nshocks)
    
    # Define our iteration function
    iterate(z_t1, ϵ_t) = T*z_t1 + R*ϵ_t
    z0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, nperiods)
        shocks[i, 1:histperiods] = histshocks[i, :]

        # Iterate state space forward
        states[:, 1, i] = iterate(z0, shocks[:, 1])
        for t in 2:nperiods
            states[:, t, i] = iterate(states[:, t-1, i], shocks[:, t])
        end

        # Apply observation and pseudo-observation equations
        observables[:, :, i]        = D        .+ Z        * states[:, :, i]
        pseudo_observables[:, :, i] = D_pseudo .+ Z_pseudo * states[:, :, i]
    end
    
    # Return shock decompositions
    return states, observables, pseudo_observables
end
