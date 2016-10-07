"""
```
smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, syses::Vector{System},
    kals::Vector{Kalman}; cond_type::Symbol = :none)

smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System}, kals::Vector{Kalman})
```

Computes and returns the smoothed values of states for every parameter draw.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses::Vector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `kals::Vector{Kalman}`: a vector of `Kalman` objects containing the results of
  the Kalman filter for each draw

### Outputs

- `alpha_hats`: a vector of smoothed states (`alpha_hat`s) returned from the
  smoother specified by `smoother_flag(m)`, one for each system in `syses`
- `eta_hats`: a vector of smoothed shocks (`eta_hat`s) returned from the
  smoother, one for each system in `syses`
"""
function smooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::DArray{System{S}, 1, Vector{System{S}}},
    kals::DArray{Kalman{S}, 1}; cond_type::Symbol = :none,
    procs::Vector{Int} = [myid()])

    data = df_to_matrix(m, df; cond_type = cond_type)
    smooth(m, data, syses, kals; procs = procs)
end

function smooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::DArray{System{S}, 1, Vector{System{S}}},
    kals::DArray{Kalman{S}, 1}; procs::Vector{Int} = [myid()])

    # numbers of useful things
    nprocs = length(procs)
    ndraws = length(syses)
    @assert length(kals) == ndraws

    nstates = n_states_augmented(m)
    nshocks = n_shocks_exogenous(m)
    nperiods = size(data, 2) - n_presample_periods(m)

    states_range = 1:nstates
    shocks_range = (nstates + 1):(nstates + nshocks)

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), procs, [nprocs])
    datas  = dfill(data, (ndraws,), procs, [nprocs])

    # Construct distributed array of smoothed states and shocks
    out = DArray((ndraws, nstates + nshocks, nperiods), procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = length(draw_inds)

        for i in draw_inds
            states, shocks = smooth(models[i], datas[i], syses[i], kals[i])

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = states
            localpart[i_local, shocks_range, :] = shocks
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:nperiods])

    return states, shocks
end

function smooth{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  sys::System{S},
                                  kal::Kalman{S})

    alpha_hat, eta_hat = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, sys, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, sys, kal[:z0], kal[:vz0])
    end

    return alpha_hat, eta_hat
end

