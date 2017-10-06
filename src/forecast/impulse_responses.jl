"""
```
impulse_responses(m, system)

impulse_responses(system, horizon)
```

Compute impulse responses for a single draw.

### Inputs

- `m::AbstractModel`: model object
- `system::System{S}`: state-space system matrices
- `horizon::Int`: number of periods ahead to forecast

where `S<:AbstractFloat`

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state impulse response functions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable impulse response functions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable impulse response functions
"""
function impulse_responses{S<:AbstractFloat}(m::AbstractModel, system::System{S})
    horizon = impulse_response_horizons(m)
    impulse_responses(system, horizon)
end

function impulse_responses{S<:AbstractFloat}(system::System{S}, horizon::Int)

    # Unpack system
    T, R = system[:TTT], system[:RRR]
    Q, Z = system[:QQ], system[:ZZ]

    Z_pseudo = system[:ZZ_pseudo]

    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 1)
    nobs         = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)


    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = zeros(S, npseudo, horizon, nshocks)

    # Define iterate function, matrix of shocks
    iterate(z_t1, ϵ_t) = T*z_t1 + R*ϵ_t
    z0 = zeros(S, nstates)
    impact = -diagm(sqrt.(diag(Q))) # a negative 1 s.d. shock

    for i = 1:nshocks
        # Iterate state space forward
        states[:, 1, i] = iterate(z0, impact[:, i])
        for t in 2:horizon
            states[:, t, i] = iterate(states[:, t-1, i], zeros(nshocks))
        end

        # Apply measurement and pseudo-measurement equations
        obs[:, :, i] = Z * states[:, :, i]
        pseudo[:, :, i] = Z_pseudo * states[:, :, i]
    end

    return states, obs, pseudo
end