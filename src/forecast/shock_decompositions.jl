"""
```
shock_decompositions{S<:AbstractFloat}(m::AbstractModel,
    systems::Vector{System{S}}, histshocks::Vector{Matrix{S}})
```

Computes shock decompositions for all draws, given a model object, system matrices, and
historical smoothed shocks.

### Inputs

- `m`: model object
- `systems::Vector{System}`: vector of length `ndraws`, whose elements are
  `System` objects specifying state-space system matrices for each draw
- `histshocks`::Vector{Matrix{S}}`: vector of length `ndraws`, whose elements
  are the `nshocks` x `hist_periods` matrices of historical smoothed shocks

### Outputs

-`states`: vector of length `ndraws`, whose elements are the `nstates` x
 `nperiods` x `nshocks` matrices of state shock decompositions
-`obs`: vector of length `ndraws`, whose elements are the `nobs` x
 `nperiods` x `nshocks` matrices of observable shock decompositions
-`pseudo_observables`: vector of length `ndraws`, whose elements are the
 `npseudo` x `nperiods` x `nshocks` matrices of pseudo-observable shock
 decompositions

where `nperiods = hist_periods + forecast_horizon`.
"""
function shock_decompositions{S<:AbstractFloat}(m::AbstractModel,
    systems::DVector{System{S}}, histshocks::DArray{S, 3};
    procs::Vector{Int} = [myid()])

    # Numbers of useful things
    ndraws = length(systems)
    nprocs = length(procs)
    horizon = forecast_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)

    # Determine periods for which to return shock decompositions
    start_ind = if !isnull(shockdec_startdate(m))
        DSGE.subtract_quarters(get(shockdec_startdate(m)), date_prezlb_start(m)) + 1
    else
        1
    end

    end_ind = if !isnull(shockdec_enddate(m))
        DSGE.subtract_quarters(get(shockdec_enddate(m)), date_prezlb_start(m)) + 1
    else
        size(histshocks, 3) + horizon
    end

    nperiods = end_ind - start_ind + 1

    # Construct distributed array of shock decompositions
    out = DArray((ndraws, nstates + nobs + npseudo, nperiods, nshocks), procs, [nprocs, 1, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds
            states, obs, pseudo = compute_shock_decompositions(systems[i], horizon,
                convert(Array, slice(histshocks, i, :, :)), start_ind, end_ind)

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :, :] = states
            localpart[i_local, obs_range,    :, :] = obs
            localpart[i_local, pseudo_range, :, :] = pseudo
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods, 1:nshocks])
    obs    = convert(DArray, out[1:ndraws, obs_range,    1:nperiods, 1:nshocks])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:nperiods, 1:nshocks])

    return states, obs, pseudo
end

"""
```
compute_shock_decompositions{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S},
    C::Vector{S}, Z::Matrix{S}, D::Vector{S}, forecast_horizons::Int,
    histshocks::Matrix{S}, start_index::Nullable{Int},
    end_index::Nullable{Int}, Z_pseudo::Matrix{S}=Matrix{S}(),
    D_pseudo::Vector{S}=Vector{S}())
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
        Array{S, 3}()
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