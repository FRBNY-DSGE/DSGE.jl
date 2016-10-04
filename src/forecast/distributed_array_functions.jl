using DataFrames, JLD

function prepare_systems_dist(m::AbstractModel, input_type::Symbol,
    params::Matrix{Float64}, TTT::Array{Float64,3}, RRR::Array{Float64,3},
    CCC::Array{Float64,3}; my_procs::Vector{Int} = [myid()])

    # Setup
    n_sim = size(params,1)
    jstep = DSGE.get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        systems = dfill(compute_system(m), (1,), my_procs)
    elseif input_type in [:full]
        empty = isempty(CCC)
        nprocs = length(my_procs);
        systems = DArray((n_sim_forecast,), my_procs, [nprocs]) do I
            draw_inds = first(I)
            ndraws_local = Int(n_sim_forecast / nprocs)
            localpart = Vector{System{Float64}}(ndraws_local)

            for i in draw_inds
                j = i * jstep
                i_local = mod(i-1, ndraws_local) + 1

                # Prepare transition eq
                TTT_j = squeeze(TTT[j, :, :], 1)
                RRR_j = squeeze(RRR[j, :, :], 1)

                if empty
                    trans_j = Transition(TTT_j, RRR_j)
                else
                    CCC_j = squeeze(CCC[j, :, :], 1)
                    trans_j = Transition(TTT_j, RRR_j, CCC_j)
                end

                # Prepare measurement eq
                params_j = vec(params[j,:])
                update!(m, params_j)
                meas_j   = measurement(m, trans_j; shocks = true)

                # Prepare system
                localpart[i_local] = System(trans_j, meas_j)
            end
            return localpart
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return systems
end


function prepare_states_dist(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
    systems::DArray{System{Float64}, 1, Vector{System{Float64}}}, params::Matrix{Float64}, df::DataFrame,
    zend::Matrix{Float64}; my_procs::Vector{Int} = [myid()])

    # Setup
    n_sim_forecast = length(systems)
    n_sim = size(params, 1)
    jstep = convert(Int, n_sim/n_sim_forecast)

    # If we just have one draw of parameters in mode, mean, or init case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        kal = DSGE.filter(m, df, systems[1]; cond_type = cond_type, allout = true)
        # `kals` is a vector of length 1
        states = dfill(kal[:filt][:, end], (1,), my_procs);

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        if cond_type in [:none]
            nprocs = length(my_procs)
            states = DArray((n_sim_forecast,), my_procs, [nprocs]) do I
                draw_inds = first(I)
                ndraws_local = Int(n_sim_forecast / nprocs)
                localpart = Vector{Vector{Float64}}(ndraws_local)

                for i in draw_inds
                    j = i * jstep
                    i_local = mod(i-1, ndraws_local) + 1

                    localpart[i_local] = vec(zend[j, :])
                end
                return localpart
            end
        elseif cond_type in [:semi, :full]
            # We will need to re-run the entire filter/smoother so we can't do anything
            # here. The reason is that while we have $s_{T|T}$ we don't have $P_{T|T}$ and
            # thus can't "restart" the Kalman filter for the conditional data period.
            states = dfill(Vector{Float64}(), (0,), my_procs)
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return states
end


function prepare_forecast_inputs_dist(m::AbstractModel, df::DataFrame;
                      input_type::Symbol  = :mode,
                      cond_type::Symbol   = :none,
                      my_procs::Vector{Int} = [myid()])
    # Some variables
    n_states = n_states_augmented(m)

    # Set up infiles
    params, TTT, RRR, CCC, zend = DSGE.load_draws(m, input_type)

    n_sim = size(params,1)
    jstep = DSGE.get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    systems = prepare_systems_dist(m, input_type, params, TTT, RRR, CCC; my_procs = my_procs)

    # Populate states vector
    states = prepare_states_dist(m, input_type, cond_type, systems, params, df, zend; my_procs = my_procs)

    return systems, states
end


function filterandsmooth_dist{T<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::DArray{System{T}, 1, Vector{System{T}}},
    z0::Vector{T} = Vector{T}(), vz0::Matrix{T} = Matrix{T}();
    cond_type::Symbol = :none, lead::Int = 0,
    my_procs::Vector{Int} = [myid()])

    data = df_to_matrix(m, df; cond_type = cond_type)
    filterandsmooth_dist(m, data, syses, z0, vz0; lead = lead, my_procs = my_procs)
end


function filterandsmooth_dist{T<:AbstractFloat}(m::AbstractModel, data::Matrix{T},
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
    zend   = DArray((ndraws,), my_procs, [nprocs]) do I
        Vector{T}[convert(Array, slice(out, i, zend_range, 1:nstates)) for i in first(I)]
    end

    # Index out SubArray for each smoothed type
    return states, shocks, pseudo, zend
end


function forecast_dist{T<:AbstractFloat}(m::AbstractModel,
    syses::DArray{System{T}, 1}, initial_state_draws::DArray{Vector{T}, 1};
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


function forecast_one_dist(m::AbstractModel{Float64}, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    output_vars::Vector{Symbol} = [], verbose::Symbol = :low,
    my_procs::Vector{Int} = [myid()])

    ### 1. Setup

    # Prepare forecast inputs
    systems, states = prepare_forecast_inputs_dist(m, df; input_type = input_type,
        cond_type = cond_type, my_procs = my_procs)

    nprocs = length(my_procs)
    ndraws = length(systems)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, DArray{Float64}}()
    forecast_output_files = DSGE.get_output_files(m, input_type, output_vars, cond_type)
    output_dir = rawpath(m, "forecast")
    # if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecasting input_type = $input_type, cond_type = $cond_type...")
        println("Forecast outputs will be saved in $output_dir")
    # end

    # Inline definition s.t. the dicts forecast_output and forecast_output_files are accessible
    function write_forecast_outputs(vars::Vector{Symbol})
        for var in vars
            file = forecast_output_files[var]
            # jldopen(file, "w") do f
            #     out = convert(Array, forecast_output[var])
            #     write(f, string(var), out)
            # end
            write_darray(file, forecast_output[var])
            # if VERBOSITY[verbose] >= VERBOSITY[:high]
                println(" * Wrote $(basename(file))")
            # end
        end
    end


    ### 2. Smoothed Histories

    # Must re-run filter/smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    filterandsmooth_vars = vcat(hist_vars, shockdec_vars)

    if !isempty(intersect(output_vars, filterandsmooth_vars)) || cond_type in [:semi, :full]

        @time histstates, histshocks, histpseudo, zends =
            filterandsmooth_dist(m, df, systems; cond_type = cond_type, my_procs = my_procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

            forecast_output[:histstates] = convert(DArray, histstates[1:end, 1:end, 1:T])
            forecast_output[:histshocks] = convert(DArray, histshocks[1:end, 1:end, 1:T])
            forecast_output[:histpseudo] = convert(DArray, histpseudo[1:end, 1:end, 1:T])
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            forecast_output[:histpseudo] = histpseudo
        end

        write_forecast_outputs(hist_vars)
    end


    ### 3. Forecasts

    # For conditional data, use the end of the hist states as the initial state
    # vector for the forecast
    if cond_type in [:semi, :full]
        states = zends
    end

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]

    if !isempty(intersect(output_vars, forecast_vars))
        @time forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast_dist(m, systems, states; my_procs = my_procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))
            histobs_cond = df_to_matrix(m, df; cond_type = cond_type)[:, index_prezlb_start(m)+T:end]
            histobs_cond = reshape(histobs_cond, (1, size(histobs_cond)...))

            nperiods = (size(histstates, 3) - T) + size(forecaststates, 3)
            nobs = n_observables(m)

            function cat_conditional(histvars::DArray{Float64, 3}, forecastvars::DArray{Float64, 3})
                nvars = size(histvars, 2)
                return DArray((ndraws, nvars, nperiods), my_procs, [nprocs, 1, 1]) do I
                    draw_inds = first(I)
                    hist_cond = convert(Array, histvars[draw_inds, :, T+1:end])
                    forecast  = convert(Array, forecastvars[draw_inds, :, :])
                    return cat(3, hist_cond, forecast)
                end
            end

            forecast_output[:forecaststates] = cat_conditional(histstates, forecaststates)
            forecast_output[:forecastshocks] = cat_conditional(histshocks, forecastshocks)
            forecast_output[:forecastpseudo] = cat_conditional(histpseudo, forecastpseudo)
            forecast_output[:forecastobs]    =
                DArray((ndraws, nobs, nperiods), my_procs, [nprocs, 1, 1]) do I
                    draw_inds = first(I)
                    ndraws_local = length(draw_inds)
                    hist_cond = repeat(histobs_cond, outer = [ndraws_local, 1, 1])
                    forecast  = convert(Array, forecastobs[draw_inds, :, :])
                    return cat(3, hist_cond, forecast)
                end
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
            forecast_output[:forecastpseudo] = forecastpseudo
            forecast_output[:forecastobs]    = forecastobs
        end

        write_forecast_outputs(forecast_vars)
    end


    ### 4. Shock Decompositions

    if !isempty(intersect(output_vars, shockdec_vars))
        @time shockdecstates, shockdecobs, shockdecpseudo =
            shock_decompositions_dist(m, systems, histshocks; my_procs = my_procs)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecpseudo] = shockdecpseudo
        forecast_output[:shockdecobs]    = shockdecobs

        write_forecast_outputs(shockdec_vars)
    end

    # Return only saved elements of dict
    filter!((k, v) -> k in output_vars, forecast_output)
    return forecast_output
end


function write_darray{T<:AbstractFloat}(filepath::AbstractString, darr::DArray{T})
    function write_localpart(pid::Int)
        jldopen(filepath, "r+") do file
            write(file, "inds$pid", collect(localindexes(darr)))
            write(file, "arr$pid", localpart(darr))
        end
    end

    jldopen(filepath, "w") do file
        write(file, "dims", darr.dims)
        write(file, "pids", collect(darr.pids))
    end

    for pid in darr.pids
        remotecall_wait(pid, write_localpart, pid)
        sleep(0.001)
    end
end

function read_darray(filepath::AbstractString)
    file = jldopen(filepath, "r")
    dims = read(file, "dims")
    pids = read(file, "pids")

    out = zeros(dims...)
    for pid in pids
        inds = read(file, "inds$pid")
        out[inds...] = read(file, "arr$pid")
    end
    close(file)
    return out
end