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
- `flip_shocks::Bool`: Whether to compute IRFs in response to a positive shock (by default the shock magnitude is a negative 1 std. shock)

where `S<:AbstractFloat`

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state impulse response functions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable impulse response functions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable impulse response functions
"""
function impulse_responses(m::AbstractRepModel, system::System{S};
                           flip_shocks::Bool = false) where {S<:AbstractFloat}
    horizon = impulse_response_horizons(m)
    states, obs, pseudo = impulse_responses(system, horizon, flip_shocks = flip_shocks)
    return states, obs, pseudo
end

function impulse_responses(system::System{S}, horizon::Int;
                           flip_shocks::Bool = false) where {S<:AbstractFloat}
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
        if flip_shocks
            shocks[i, 1] = sqrt(system[:QQ][i, i]) # a positive 1 s.d. shock
        else
            shocks[i, 1] = -sqrt(system[:QQ][i, i]) # a negative 1 s.d. shock
        end
        # Iterate state space forward
        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(system, s_0, shocks)
    end

    return states, obs, pseudo
end

# Method for specifying the subset of shocks, and the size of each shock
function impulse_responses(m::AbstractModel, system::System{S},
                           horizon::Int, shock_names::Vector{Symbol},
                           shock_values::Vector{Float64}) where S<:AbstractFloat

    # Must provide a name and value for each shock
    @assert length(shock_names) == length(shock_values)

    # Setup
    exo          = m.exogenous_shocks
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

    for (i, shock) in enumerate(shock_names)
        # Isolate single shock
        shocks = zeros(S, nshocks, horizon)
        shocks[exo[shock], 1] = shock_values[i]

        # Iterate state space forward
        states[:, :, exo[shock]], obs[:, :, exo[shock]], pseudo[:, :, exo[shock]], _ = forecast(system, s_0, shocks)
    end

    return states, obs, pseudo
end

# Method for picking a specific shock and the size of the desired shift in the state or
# observed variable using that shock. (Omitting the feature to do pseudo-observables
# now, for simplicity)
function impulse_responses(m::AbstractModel, system::System{S},
                           horizon::Int, shock_name::Symbol,
                           var_name::Symbol, var_value::Float64) where S<:AbstractFloat
    # Setup
    var_names, var_class =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states
    elseif var_name in keys(m.observables)
        m.observables, :obs
    end
    exo          = m.exogenous_shocks
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

    # Isolate single shock
    shocks = zeros(S, nshocks, horizon)
    if var_class == :states
        shocks[exo[shock_name], 1] = obtain_shock_from_desired_state_value(var_value,
                                                                      var_names[var_name],
                                                                      exo[shock_name],
                                                                      system[:RRR])
    else # == :obs
        shocks[exo[shock_name], 1] = obtain_shock_from_desired_obs_value(var_value,
                                                                    var_names[var_name],
                                                                    exo[shock_name],
                                                                    system[:ZZ],
                                                                    system[:RRR])
    end

    # Iterate state space forward
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)

    return states, obs, pseudo
end

# Given a desired value, x, for some state s^i_1, back out the necessary value of
# ϵ^j_1 needed s.t. s^i_1 = x.
# E.g. if I want the impulse response of output (s^i_1) to be 50 basis points
# based on a shock to productivity (ϵ^j_1), how large does that shock
# to productivity need to be?
# The answer is ϵ^j_1 = x/R_{i,j}
function obtain_shock_from_desired_state_value(state_value::Float64, state_ind::Int,
                                               shock_ind::Int, RRR::Matrix{Float64})
    return state_value/RRR[state_ind, shock_ind]
end

# Similar principle to the function for states
# The answer is ϵ^j_1 = x/{Σ_{k=1}^n Z_{i, k} R_{k, j}},
# where n is the total number of states, i is the index of the observable of interest
# y^i_1, and j is the index of the shock of interest, ϵ^j_1
function obtain_shock_from_desired_obs_value(obs_value::Float64, obs_ind::Int,
                                             shock_ind::Int, ZZ::Matrix{Float64},
                                             RRR::Matrix{Float64})
    return obs_value/dot(ZZ[obs_ind, :], RRR[:, shock_ind])
end
