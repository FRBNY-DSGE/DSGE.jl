function filterandsmooth_dist{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
                                                syses::DArray{System{S}, 1},
                                                z0::Vector{S} = Vector{S}(),
                                                vz0::Matrix{S} = Matrix{S}();
                                                lead::Int = 0,
                                                my_procs::Vector{Int} = [myid()])

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
        ndraws_local = Int(ndraws / nprocs)

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
    zend   = convert(DArray, out[1:ndraws, zend_range,   1:nstates])

    # Index out SubArray for each smoothed type
    return states, shocks, pseudo, zend
end


function forecast_dist{T<:AbstractFloat}(m::AbstractModel, syses::DArray{System{T}, 1},
                                         initial_state_draws::DArray{Vector{T}, 1};
                                         shock_distributions::Union{Distribution,Matrix{T}} = Matrix{T}(),
                                         my_procs::Vector{Int} = [myid()])


    # Numbers of useful things
    ndraws = length(syses)
    nprocs = length(my_procs)
    horizon = forecast_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = 12
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)
    shocks_range = (nstates + nobs + npseudo + 1):(nstates + nobs + npseudo + nshocks)

    # set up distribution of shocks if not specified
    # For now, we construct a giant vector of distirbutions of shocks and pass
    # each to compute_forecast.
    #
    # TODO: refactor so that compute_forecast
    # creates its own DegenerateMvNormal based on passing the QQ
    # matrix (which has already been computed/is taking up space)
    # rather than having to copy each Distribution across nodes. This will also be much more
    # space-efficient when forecast_kill_shocks is true.

    shock_distributions = if isempty(shock_distributions)
        if forecast_kill_shocks(m)
            dfill(zeros(nshocks, horizon), (ndraws,), my_procs, [nprocs])
        else
            # use t-distributed shocks
            if forecast_tdist_shocks(m)
                dfill(Distributions.TDist(forecast_tdist_df_val(m)), (ndraws,), my_procs, [nprocs])
            # use normally distributed shocks
            else
                DArray(I -> [DSGE.DegenerateMvNormal(zeros(nshocks), sqrt(s[:QQ])) for s in syses[I...]],
                       (ndraws,), my_procs, [nprocs])
            end
        end
    end

    # Construct distributed array of forecast outputs
    out = DArray((ndraws, nstates + nobs + npseudo + nshocks, horizon), my_procs, [nprocs, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds
            dict = DSGE.compute_forecast(syses[i], horizon, shock_distributions[i], initial_state_draws[i])

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :] = dict[:states]
            localpart[i_local, obs_range,    :] = dict[:observables]
            localpart[i_local, pseudo_range, :] = dict[:pseudo_observables]
            localpart[i_local, shocks_range, :] = dict[:shocks]
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:horizon])
    obs    = convert(DArray, out[1:ndraws, obs_range,    1:horizon])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:horizon])
    shocks = convert(DArray, out[1:ndraws, shocks_range, 1:horizon])

    return states, obs, pseudo, shocks
end


function shock_decompositions_dist{T<:AbstractFloat}(m::AbstractModel,
    syses::DArray{System{T}, 1}, histshocks::DArray{T, 3};
    my_procs::Vector{Int} = [myid()])

    # Numbers of useful things
    ndraws = length(syses)
    nprocs = length(my_procs)
    horizon = forecast_horizons(m)

    nstates = n_states_augmented(m)
    nobs    = n_observables(m)
    npseudo = 12
    nshocks = n_shocks_exogenous(m)

    states_range = 1:nstates
    obs_range    = (nstates + 1):(nstates + nobs)
    pseudo_range = (nstates + nobs + 1):(nstates + nobs + npseudo)

    # Determine periods for which to return shock decompositions
    start_ind = if !isnull(shockdec_startdate(m))
        DSGE.subtract_quarters(get(shockdec_startdate(m)), date_prezlb_start(m)) + 1
    else
        1
    end

    end_ind = if !isnull(shockdec_enddate(m))
        DSGE.subtract_quarters(get(shockdec_enddate(m)), date_prezlb_start(m)) + 1
    else
        size(histshocks, 3) + horizon
    end

    nperiods = end_ind - start_ind + 1

    # Construct distributed array of shock decompositions
    out = DArray((ndraws, nstates + nobs + npseudo, nperiods, nshocks), my_procs, [nprocs, 1, 1, 1]) do I
        localpart = zeros(map(length, I)...)
        draw_inds = first(I)
        ndraws_local = Int(ndraws / nprocs)

        for i in draw_inds
            states, obs, pseudo = DSGE.compute_shock_decompositions(syses[i], horizon,
                convert(Array, slice(histshocks, i, :, :)), start_ind, end_ind)

            i_local = mod(i-1, ndraws_local) + 1

            localpart[i_local, states_range, :, :] = states
            localpart[i_local, obs_range,    :, :] = obs
            localpart[i_local, pseudo_range, :, :] = pseudo
        end
        return localpart
    end

    # Convert SubArrays to DArrays and return
    states = convert(DArray, out[1:ndraws, states_range, 1:nperiods, 1:nshocks])
    obs    = convert(DArray, out[1:ndraws, obs_range,    1:nperiods, 1:nshocks])
    pseudo = convert(DArray, out[1:ndraws, pseudo_range, 1:nperiods, 1:nshocks])

    return states, obs, pseudo
end
