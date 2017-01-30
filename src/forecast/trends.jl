"""
```
trends(m, systems; procs = [myid()])
```

Computes trend (steady-state) values of states, observables, and
pseudoobservables, given a model object and system matrices.

### Inputs

- `m::AbstractModel`: model object
- `systems::Vector{System{S}}`: vector of `System` objects specifying
  state-space system matrices for each draw

where `S<:AbstractFloat`.

### Outputs

- `states::Matrix{S}`: matrix of size `ndraws` x `nstates` of state
  steady-state values for each draw
- `obs::Matrix{S}`: matrix of size `ndraws` x `nobs` of observable
  steady-state for each draw
- `pseudo::Matrix{S}`: matrix of size `ndraws` x `npseudo` of
  pseudo-observable steady-state values for each draw. If
  `!forecast_pseudoobservables(m)`, `pseudo` will be empty.
"""
function trends{S<:AbstractFloat}(m::AbstractModel{S}, systems::Vector{System{S}})

    # Numbers of useful things
    ndraws  = length(systems)
    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = n_pseudoobservables(m)

    states = zeros(ndraws, nstates)
    obs    = zeros(ndraws, nobs)
    pseudo = zeros(ndraws, npseudo)

    for i = 1:ndraws
        states_i, obs_i, pseudo_i = compute_trend(systems[i])

        states[i, :] = states_i
        obs[i,    :] = obs_i
        pseudo[i, :] = pseudo_i
    end

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
