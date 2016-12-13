"""
```
deterministic_trends(m, systems, z0s; procs = [myid()])
```

Computes deterministic trend values of states, observables, and
pseudoobservables for all draws, given a model object and system
matrices. The deterministic trend for a single draw is simply the
series that would be obtained by iterating the state-space system
forward, beginning from a state vector z0 in the first historical period.

### Inputs

- `m::AbstractModel`: model object
- `systems::DVector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw
- `z0s`::DVector{Vector{Float64}}: vector of initial state vectors for each draw

where `S<:AbstractFloat`.

### Keyword Arguments

- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`

### Outputs

- `states::DArray{S, 3}`: array of size `ndraws` x `nstates` x `nperiods` of state
  steady-state values for each draw
- `obs::DArray{S, 3}`: array of size `ndraws` x `nobs` x `nperiods` of observable
  steady-state for each draw
- `pseudo::DArray{S, 3}`: array of size `ndraws` x `npseudo` x `nperiods` of
  pseudo-observable steady-state values for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.

where `nperiods` is the number of quarters between
`date_mainsample_start(m)` and `date_forecast_end(m)`, inclusive.
"""
function deterministic_trends{S<:AbstractFloat}(m::AbstractModel, systems::DVector{System{S}},
                                                z0s::DVector{Vector{S}}; procs::Vector{Int} = [myid()])


    @assert length(systems) == length(z0s)
    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)
    nprocs = length(procs)

    # Numbers of useful things
    ndraws  = length(systems)
    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods  = DSGE.subtract_quarters(date_shockdec_end(m), date_mainsample_start(m)) + 1
    n_presample = n_presample_periods(m)
    date_range = (index_shockdec_start(m) - n_presample):(index_shockdec_end(m) - n_presample)

    # Assign ranges for DArray
    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)

    # Compute deterministic trend for each draw, return DArray
    out = DArray((ndraws, nstates + nobs + npseudo, nperiods),
                       procs, [nprocs, 1, 1]) do I

        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)   # Number of draws that localpart stores data for

        for i in draw_inds
            dettrend_states, dettrend_obs, dettrend_pseudo =
                compute_deterministic_trend(systems[i], z0s[i], nperiods)

            # Assign the i-th index of systems (the draw)
            # to the i_local-th index of the localpart array
            i_local = mod(i-1, ndraws_local) + 1

            # Assign return values from compute_shock_decompositions to a slice of localpart
            localpart[i_local, states_range, :, :] = dettrend_states
            localpart[i_local, obs_range,    :, :] = dettrend_obs
            if forecast_pseudoobservables(m)
                localpart[i_local, pseudo_range, :, :] = dettrend_pseudo
            end
        end
        return localpart
    end

    # Convert SubArrays (returned when indexing `out`) to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, date_range])
    obs    = convert(DArray, out[1:ndraws, obs_range,    date_range])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, date_range])

    return states, obs, pseudo
end


function compute_deterministic_trend{S<:AbstractFloat}(system::System{S}, z0::Vector{S}, nperiods::Int)

    # construct matrix of 0 shocks for entire history and forecast horizon
    nshocks  = size(system[:RRR], 2)
    shocks   = zeros(S, nshocks, nperiods)

    # use compute_forecast to iterate state-space system forward without shocks or the ZLB procedure
    states, obs, pseudo, _ = compute_forecast(system, z0, shocks)

    states, obs, pseudo
end
