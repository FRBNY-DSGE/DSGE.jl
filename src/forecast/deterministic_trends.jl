"""
```
deterministic_trends(m, systems, z0s)
```

Computes deterministic trend values of states, observables, and
pseudoobservables for all draws, given a model object and system matrices. The
deterministic trend for a single draw is simply the series that would be
obtained by iterating the state-space system forward, beginning from a state
vector `z0` in the last presample period.

### Inputs

- `m::AbstractModel`: model object
- `systems::Vector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw
- `z0s`::Vector{Vector{S}}: vector of initial state vectors for each draw

where `S<:AbstractFloat`.

### Outputs

- `states::Array{S, 3}`: array of size `ndraws` x `nstates` x `nperiods` of state
  steady-state values for each draw
- `obs::Array{S, 3}`: array of size `ndraws` x `nobs` x `nperiods` of observable
  steady-state for each draw
- `pseudo::Array{S, 3}`: array of size `ndraws` x `npseudo` x `nperiods` of
  pseudo-observable steady-state values for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.

where `nperiods` is the number of quarters between
`date_mainsample_start(m)` and `date_forecast_end(m)`, inclusive.
"""
function deterministic_trends{S<:AbstractFloat}(m::AbstractModel{S},
    systems::Vector{System{S}}, z0s::Vector{Vector{S}})

    @assert length(systems) == length(z0s)

    # Numbers of useful things
    ndraws  = length(systems)
    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)

    # Dates: We compute the deterministic trend starting from the
    # first historical period.  However, since it is only used to
    # compute shock decompositions, we truncate and only store
    # results for periods corresponding to the shockdec period.
    nperiods   = subtract_quarters(date_shockdec_end(m), date_mainsample_start(m)) + 1
    date_range = index_shockdec_start(m):index_shockdec_end(m)

    states = zeros(ndraws, nstates, nperiods)
    obs    = zeros(ndraws, nobs,    nperiods)
    pseudo = zeros(ndraws, npseudo, nperiods)

    for i = 1:ndraws
        states_i, obs_i, pseudo_i = compute_deterministic_trend(systems[i], z0s[i], nperiods)

        states[i, :, :] = states_i
        obs[i,    :, :] = obs_i
        pseudo[i, :, :] = pseudo_i
    end

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
