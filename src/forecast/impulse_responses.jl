"""
```
impulse_responses(m, systems; procs = [myid()])
```

Computes impulse responses for all states, pseudo-observables, or observables
across all draws.

### Inputs
- `m::AbstractModel`: model object
- `systems::DVector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw
- `impulse_response_shocks`: a Dict mapping the desired shocks to their indices

### Keyword Arguments

- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`

### Outputs

- `states::DArray{S, 4}`: array of size `ndraws` x `nstates` x `horizon` x
  `nshocks` of state shock decompositions for each draw
- `obs::DArray{S, 4}`: array of size `ndraws` x `nobs` x `horizon` x `nshocks`
  of observable shock decompositions for each draw
- `pseudo::DArray{S, 4}`: array of size `ndraws` x `npseudo` x `horizon` x
  `nshocks` of pseudo-observable shock decompositions for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.

where `horizon` is the forecast horizon for the model as given by
`impulse_response_horizons(m)`
"""
function impulse_responses{S<:AbstractFloat}(m::AbstractModel,
                                             systems::DVector{System{S}, Vector{System{S}}};
                                             procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs)

    # Numbers of useful things
    ndraws = length(systems)
    nprocs = length(procs)
    horizon = impulse_response_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)

    # Construct distributed array of shock decompositions
    out = DArray((ndraws, nstates + nobs + npseudo, horizon, nshocks), procs, [nprocs, 1, 1, 1]) do I # I is a tuple of indices for given localpart

        # Compute shock decomposition for each draw
        localpart = zeros(map(length, I)...)  # Regular array. In this DArray constructor, each process has one localpart.
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)   # Number of draws that localpart stores data for

        for i in draw_inds
            states, obs, pseudo = compute_impulse_response(systems[i], horizon)

            # Assign the i-th index of systems (the draw)
            # to the i_local-th index of the localpart array
            i_local = mod(i-1, ndraws_local) + 1

            # Assign return values from compute_shock_decompositions to a slice of localpart
            localpart[i_local, states_range, :, :] = states
            localpart[i_local, obs_range,    :, :] = obs
            if forecast_pseudoobservables(m)
                localpart[i_local, pseudo_range, :, :] = pseudo
            end
        end
        return localpart
    end

    # Convert SubArrays (returned when indexing `out`) to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:horizon, 1:nshocks])
    obs    = convert(DArray, out[1:ndraws, obs_range,    1:horizon, 1:nshocks])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:horizon, 1:nshocks])

    return states, obs, pseudo
end


"""
```
compute_impulse_response(system, horizon)
```

Compute impulse responses for a single draw.

### Inputs

- `system::System{S}`: state-space system matrices. Alternatively, provide
  transition equation matrices `T`, `R`; measurement equation matrices `Z`, `D`;
  and (possibly empty) pseudo-measurement equation matrices `Z_pseudo` and
  `D_pseudo`.
- `horizon::Int`: number of periods ahead to forecast

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state shock decompositions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable shock decompositions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable shock decompositions. If the provided `Z_pseudo` and
  `D_pseudo` matrices are empty, then `pseudo` will be empty.
"""
function compute_impulse_response{S<:AbstractFloat}(system::System{S}, horizon::Int)

    # Unpack system
    T, R = system[:TTT], system[:RRR]
    Q, Z = system[:QQ], system[:ZZ]

    Z_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo]
    else
        Matrix{S}()
    end

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 1)
    nobs         = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)

    forecast_pseudo = !isempty(Z_pseudo)

    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = if forecast_pseudo
        zeros(S, npseudo, horizon, nshocks)
    else
        zeros(S, 0, 0, 0)
    end

    # Define iterate function, matrix of shocks
    iterate(z_t1, ϵ_t) = T*z_t1 + R*ϵ_t
    z0 = zeros(S, nstates)
    impact = -diagm(sqrt(diag(Q))) # a negative 1 s.d. shock

    for i = 1:nshocks
        # Iterate state space forward
        states[:, 1, i] = iterate(z0, impact[:, i])
        for t in 2:horizon
            states[:, t, i] = iterate(states[:, t-1, i], zeros(nshocks))
        end

        # Apply measurement and pseudo-measurement equations
        obs[:, :, i] = Z * states[:, :, i]
        if forecast_pseudo
            pseudo[:, :, i] = Z_pseudo * states[:, :, i]
        end
    end

    return states, obs, pseudo
end