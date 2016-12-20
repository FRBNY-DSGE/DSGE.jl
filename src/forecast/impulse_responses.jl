"""
```
impulse_responses(m, systems; procs = [myid()])
impulse_responses(m, systems, impulse_response_shocks; procs = [myid()])
```

computes impulse responses for all states, pseudo-observables, or observables
for a set of shocks across all draws. 

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
    `impulse_response_horizons(m)` and `nshocks` is the number of shocks in
    `impulse_response_shocks`.

"""
function impulse_responses{S<:AbstractFloat}(m::AbstractModel,
                                systems::DVector{System{S}, Vector{System{S}}};
                                procs::Vector{Int} = [myid()])

    impulse_response_shocks = m.exogenous_shocks

    return impulse_responses(m, systems, impulse_response_shocks; procs = procs)
end

function impulse_responses{S<:AbstractFloat}(m::AbstractModel,
                                systems::DVector{System{S}, Vector{System{S}}},
                                impulse_response_shocks::Dict{Symbol,Int64};
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
    nshocks = length(impulse_response_shocks)

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
            states, obs, pseudo = compute_impulse_response(systems[i], horizon, impulse_response_shocks)

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
compute_impulse_response(system, horizon, impulse_response_shocks)

compute_impulse_response(T, R, Z, D, Z_pseudo, D_pseudo, QQ, horizon, impulse_response_shocks)
```

compute impulse responses for a single 

### Inputs

- `system::System{S}`: state-space system matrices. Alternatively, provide
  transition equation matrices `T`, `R`; measurement equation matrices `Z`, `D`;
  and (possibly empty) pseudo-measurement equation matrices `Z_pseudo` and
  `D_pseudo`.
- `horizon::Int`: number of periods ahead to forecast
- `impulse_response_shocks`: Dictionary mapping shocks to indices.
   where `S<:AbstractFloat`.

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state shock decompositions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable shock decompositions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable shock decompositions. If the provided `Z_pseudo` and
  `D_pseudo` matrices are empty, then `pseudo` will be empty.

Only the impulse_responses corresponding to those in `impulse_response_shocks` will be calculated and returned. 
"""
function compute_impulse_response{S<:AbstractFloat}(system::System{S},
    horizon::Int, impulse_response_shocks::Dict{Symbol,Int64})

    # Unpack system
    T, R = system[:TTT], system[:RRR]
    Z, D = system[:ZZ], system[:DD]
    Q    = system[:QQ]

    Z_pseudo, D_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    compute_impulse_response(T, R, Z, D, Z_pseudo, D_pseudo,
        Q, horizon, impulse_response_shocks)
end

function compute_impulse_response{S<:AbstractFloat}(T::Matrix{S},
    R::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, Q::Matrix{S}, horizon::Int,
    impulse_response_shocks::Dict{Symbol,Int64})

    # Setup
    nshocks      = size(Q, 1)
    nstates      = size(T, 2)
    nobs         = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
 
    forecast_pseudo = !isempty(Z_pseudo) && !isempty(D_pseudo)

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
    shock_indices = collect(values(impulse_response_shocks))
    impact = -diagm(sqrt(diag(Q))) # a negative 1 s.d. shock

    for i in shock_indices               
        # Iterate state space forward
        states[:, 1, i] = iterate(z0, impact[:,i])
        for t in 2:horizon
            states[:, t, i] = iterate(states[:, t-1, i], zeros(size(Q,1)))
        end

        # Apply measurement and pseudo-measurement equations
        obs[:, :, i] = Z * states[:, :, i]
        if forecast_pseudo
            pseudo[:, :, i] = Z_pseudo * states[:, :, i]
        end
    end

    # states = states[:,:,shock_indices]
    # obs    = obs[:,:,shock_indices]
    # if forecast_pseudo
    #     pseudo = pseudo[:,:,shock_indices]
    # end
    
    return states, obs, pseudo
end


