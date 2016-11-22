"""
```
trends(m, systems; procs = [myid()])
```

Computes trend (steady-state) values of states, observables, and
pseudoobservables, given a model object and system matrices.

### Inputs

- `m::AbstractModel`: model object
- `systems::DVector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw

where `S<:AbstractFloat`.

### Keyword Arguments

- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`

### Outputs

- `states::DArray{S, 2}`: DArray of size `ndraws` x `nstates` of state
  steady-state values for each draw
- `obs::DArray{S, 2}`: DArray of size `ndraws` x `nobs` of observable
  steady-state for each draw
- `pseudo::DArray{S, 2}`: DArray of size `ndraws` x `npseudo` of
  pseudo-observable steady-state values for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
"""
function trends{S<:AbstractFloat}(m::AbstractModel, systems::DVector{System{S}};
                                  procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)
    nprocs = length(procs)

    # Numbers of useful things
    ndraws  = length(systems)
    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    # Assign ranges for DArray
    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)

    # Compute trend for each draw, return DArray
    # Return distributed array of trends for all draws
    out = DArray((ndraws, nstates + nobs + npseudo), procs, [nprocs, 1]) do I

        # Construct localpart array
        localpart = zeros(map(length, I))
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds

            # Compute trend
            states_trend, obs_trend, pseudo_trend = compute_trend(systems[i])

            # Assign the i-th index of systems (the draw)
            # to the i_local-th index of the localpart array
            i_local = mod(i-1, ndraws_local) + 1

            # Assign return values from compute_shock_decompositions to a slice of localpart
            localpart[i_local, states_range] = states_trend
            localpart[i_local, obs_range] = obs_trend
            if npseudo > 0
                localpart[i_local, pseudo_range] = pseudo_trend
            end
        end
        return localpart
    end

    # Convert SubArrays (returned when indexing `out`) to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range])
    obs    = convert(DArray, out[1:ndraws, obs_range   ])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range])

    return states, obs, pseudo
end

"""
```
compute_trend{S<:AbstractFloat}(system::System{S})
```

Compute trend (steady-state) states, observables, and
pseudoobservables. The trend is used for plotting shock
decompositions.
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
