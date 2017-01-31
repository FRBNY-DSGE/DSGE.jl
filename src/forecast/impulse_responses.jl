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

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state impulse response functions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable impulse response functions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable impulse response functions. If the pseudo-measurement equation
  matrices in `system` are empty, then `pseudo` will be empty.
"""
function impulse_responses{S<:AbstractFloat}(m::AbstractModel, system::System{S})
    horizon = impulse_response_horizons(m)
    impulse_responses(system, horizon)
end

function impulse_responses{S<:AbstractFloat}(system::System{S}, horizon::Int)

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