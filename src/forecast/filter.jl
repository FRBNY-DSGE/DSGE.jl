"""
```
filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
    include_presample::Bool = true)

filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)
```

Computes and returns the filtered values of states for every state-space system in `syses`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `z0`: an optional `Nz` x 1 initial state vector
- `vz0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector
- `allout`: an optional keyword argument indicating whether we want optional
  output variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned vector of `Kalman` objects

### Outputs

`filter` returns a vector of `Kalman` objects, which each contain the following fields:

- `logl`: value of the average log likelihood function of the SSM under assumption that
  observation noise ϵ(t) is normally distributed
- `pred`: a `Nz` x `T+lead` matrix containing one-step predicted state vectors.
- `vpred`: a `Nz` x `Nz` x `T+lead` matrix containing mean square errors of predicted
  state vectors.
- `filt`: an optional `Nz` x `T` matrix containing filtered state vectors.
- `vfilt`: an optional `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
  state vectors.
"""
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::DArray{System{S}, 1, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
    include_presample::Bool = true, my_procs::Vector{Int} = [myid()])

    # Convert the DataFrame to a data matrix without altering the original dataframe
    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, syses, z0, vz0; lead = lead, allout = allout,
           include_presample = include_presample, my_procs = my_procs)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::DArray{System{S}, 1, Vector{System{S}}},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
    lead::Int = 0, allout::Bool = false, include_presample::Bool = true,
    my_procs::Vector{Int} = [myid()])

    # Numbers of useful things
    ndraws = length(syses)
    nprocs = length(my_procs)

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), my_procs, [nprocs])
    datas  = dfill(data, (ndraws,), my_procs, [nprocs])
    z0s    = dfill(z0,   (ndraws,), my_procs, [nprocs])
    vz0s   = dfill(vz0,  (ndraws,), my_procs, [nprocs])

    # Construct distributed array of Kalman objects
    kals = DArray((ndraws,), my_procs, [nprocs]) do I
        draw_inds = first(I)
        ndraws_local = length(draw_inds)
        localpart = Vector{Kalman{S}}(ndraws_local)

        for i in draw_inds
            i_local = mod(i-1, ndraws_local) + 1
            localpart[i_local] = filter(models[i], datas[i], syses[i], z0s[i], vz0s[i];
                                        allout = allout, include_presample = include_presample)
        end
        return localpart
    end
end

function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, sys::System{S},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  cond_type::Symbol = :none, lead::Int = 0,
                                  allout::Bool = false, include_presample::Bool = true)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, sys, z0, vz0; lead = lead, allout = allout, include_presample = include_presample)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System{S},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  lead::Int = 0, allout::Bool = false, include_presample::Bool = true)

    # pull out the elements of sys
    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    DD     = sys[:DD]
    VVall  = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, TTT, CCC, ZZ, DD, VVall, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    end

    return kal
end

"""
```
filterandsmooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)

filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)
```

Computes and returns the smoothed states, shocks, and pseudo-observables, as
well as the Kalman filter outputs, for every state-space system in `syses`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `z0`: an optional `Nz` x 1 initial state vector
- `vz0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector

### Outputs

- `states`: 3-dimensional array of size `nstates` x `hist_periods` x `ndraws`
  consisting of smoothed states for each draw
- `shocks`: 3-dimensional array of size `nshocks` x `hist_periods` x `ndraws`
  consisting of smoothed shocks for each draw
- `pseudo`: 3-dimensional array of size `npseudo` x `hist_periods` x `ndraws`
  consisting of pseudo-observables computed from the smoothed states for each
  draw
- `kals`: vector of Kalman objects, of length `ndraws`

where `states` and `shocks` are returned from the smoother specified by
`smoother_flag(m)`.
"""
function filterandsmooth{T<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::DArray{System{T}, 1, Vector{System{T}}},
    z0::Vector{T} = Vector{T}(), vz0::Matrix{T} = Matrix{T}();
    cond_type::Symbol = :none, lead::Int = 0,
    my_procs::Vector{Int} = [myid()])

    data = df_to_matrix(m, df; cond_type = cond_type)
    filterandsmooth(m, data, syses, z0, vz0; lead = lead, my_procs = my_procs)
end


function filterandsmooth{T<:AbstractFloat}(m::AbstractModel, data::Matrix{T},
    syses::DArray{System{T}, 1, Vector{System{T}}},
    z0::Vector{T} = Vector{T}(), vz0::Matrix{T} = Matrix{T}();
    lead::Int = 0, my_procs::Vector{Int} = [myid()])

    # Numbers of useful things
    ndraws = length(syses)
    nprocs = length(my_procs)
    nperiods = size(data, 2) - n_presample_periods(m)

    nstates = n_states_augmented(m)
    npseudo = 12
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    shocks_range = (nstates + 1):(nstates + nshocks)
    pseudo_range = (nstates + nshocks + 1):(nstates + nshocks + npseudo)
    zend_range   = nstates + nshocks + npseudo + 1

    # Broadcast models and data matrices
    models = dfill(m,    (ndraws,), my_procs, [nprocs])
    datas  = dfill(data, (ndraws,), my_procs, [nprocs])
    z0s    = dfill(z0,   (ndraws,), my_procs, [nprocs])
    vz0s   = dfill(vz0,  (ndraws,), my_procs, [nprocs])

    # Construct distributed array of smoothed states, shocks, and pseudo-observables
    out = DArray((ndraws, nstates + nshocks + npseudo + 1, nperiods), my_procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = length(draw_inds)

        for i in draw_inds
            states, shocks, pseudo, zend = filterandsmooth(models[i], datas[i], syses[i], z0s[i], vz0s[i])

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = states
            localpart[i_local, shocks_range, :] = shocks
            localpart[i_local, pseudo_range, :] = pseudo
            localpart[i_local, zend_range,   1:nstates] = zend
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:nperiods])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:nperiods])
    zend   = DArray((ndraws,), my_procs, [nprocs]) do I
        Vector{T}[convert(Array, slice(out, i, zend_range, 1:nstates)) for i in first(I)]
    end

    # Index out SubArray for each smoothed type
    return states, shocks, pseudo, zend
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System{S},
                                           z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                           lead::Int = 0)

    TTT   = sys[:TTT]
    RRR   = sys[:RRR]
    CCC   = sys[:CCC]
    QQ    = sys[:QQ]
    ZZ    = sys[:ZZ]
    DD    = sys[:DD]
    VVall = sys[:VVall]

    filterandsmooth(m, data, TTT, RRR, CCC, QQ, ZZ, DD, VVall, z0, vz0; lead = lead)
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
                                           TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
                                           QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S}, VVall::Matrix{S},
                                           z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                           lead::Int = 0)
    ## 1. Filter

    # Call the appropriate version of the Kalman filter
    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0; lead =
            lead, allout = true, include_presample = true)
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, TTT, CCC, ZZ, DD, VVall, z0, vz0;
            lead = lead, allout = true, include_presample = true)
    end

    ## 2. Smooth

    states, shocks = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD,
                        kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD,
                                kal[:z0], kal[:vz0])
    end

    ## 3. Map smoothed states to pseudo-observables

    # for now, we are ignoring pseudo-observables so these can be empty
    Z_pseudo = zeros(S, 12, n_states_augmented(m))
    D_pseudo = zeros(S, 12)

    pseudo = D_pseudo .+ Z_pseudo * states

    return states, shocks, pseudo, kal[:zend]
end
