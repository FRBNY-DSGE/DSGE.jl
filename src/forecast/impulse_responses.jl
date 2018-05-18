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
    # Setup
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = zeros(S, npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    s_0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, horizon)
        shocks[i, 1] = -sqrt(system[:QQ][i, i]) # a negative 1 s.d. shock

        # Iterate state space forward
        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(system, s_0, shocks)
    end

    return states, obs, pseudo
end