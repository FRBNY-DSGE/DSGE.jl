function decompose_forecast_new(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                                params_new::Vector{Float64}, params_old::Vector{Float64},
                                cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol};
                                individual_shocks::Bool = false,
                                check::Bool = false, atol::Float64 = 1e-8) where M<:AbstractModel
    # Compute numbers of periods
    T0, T, k, T1_new, T1_old, _ = DSGE.decomposition_periods(m_new, m_old, cond_new, cond_old)
    @assert size(df_new, 1) == T0 + T + T1_new
    @assert size(df_old, 1) == T0 + T - k + T1_old

    # Forecast
    keep_startdate = date_forecast_start(m_new) # T+1
    keep_enddate   = date_forecast_end(m_old)   # T+H
    H = DSGE.subtract_quarters(keep_enddate, keep_startdate) + 1
    shockdec_splitdate = date_mainsample_end(m_old) # T-k

    f(m::AbstractModel, df::DataFrame, params::Vector{Float64}, cond_type::Symbol; kwargs...) =
        decomposition_forecast(m, df, params, cond_type, keep_startdate, keep_enddate,
                               shockdec_splitdate; kwargs..., check = check)

    out1 = f(m_new, df_new, params_new, cond_new, outputs = [:shockdec])            # new data, new params
    out2 = f(m_old, df_old, params_new, cond_old, outputs = [:shockdec, :forecast]) # old data, new params
    out3 = f(m_old, df_old, params_old, cond_old, outputs = [:forecast])            # old data, old params

    # Initialize output dictionary
    decomp = Dict{Symbol, Array{Float64}}()

    # Decomposition
    for class in classes
        # All elements of out are of size Ny x Nh, where the second dimension
        # ranges from t = T+1:T+H
        forecastvar = Symbol(:forecast, class) # Z s_{T+h} + D
        trendvar    = Symbol(:trend,    class) # Z D
        dettrendvar = Symbol(:dettrend, class) # Z T^{T+h} s_0 + D
        revisionvar = Symbol(:revision, class) # Z \sum_{t=1}^{T-k} T^{T+h-t} R ϵ_t + D
        newsvar     = Symbol(:news,     class) # Z \sum_{t=T-k+1}^{T+h} T^{T+h-t} R ϵ_t + D

        # Data revision
        data_comp = (out1[dettrendvar] - out2[dettrendvar]) + (out1[revisionvar] - out2[revisionvar])
        decomp[Symbol(:decompdata, class)] = data_comp

        # News
        news_comp = out1[newsvar] - out2[newsvar]
        decomp[Symbol(:decompnews, class)] = news_comp

        # Parameter re-estimation
        para_comp = out2[forecastvar] - out3[forecastvar]
        decomp[Symbol(:decomppara, class)] = para_comp

        # Check decomposition
        if check
            exp_diff = out1[forecastvar] - out3[forecastvar]
            act_diff = para_comp + data_comp + news_comp
            @assert act_diff ≈ exp_diff
        end
    end

    return decomp, out1, out2, out3
end

function decomposition_forecast(m::AbstractModel, df::DataFrame, params::Vector{Float64}, cond_type::Symbol,
                                keep_startdate::Date, keep_enddate::Date, shockdec_splitdate::Date;
                                outputs::Vector{Symbol} = [:forecast, :shockdec], check::Bool = false)
    # Compute state space
    DSGE.update!(m, params)
    system = compute_system(m)

    # Initialize output dictionary
    out = Dict{Symbol, Array{Float64}}()

    # Smooth
    histstates, histshocks, histpseudo, s_0 = smooth(m, df, system, cond_type = cond_type, draw_states = false)

    # Forecast
    if :forecast in outputs || check
        s_T = histstates[:, end]
        _, forecastobs, forecastpseudo, _ =
            forecast(m, system, s_T, cond_type = cond_type, enforce_zlb = false, draw_shocks = false)

        T = n_mainsample_periods(m)
        if cond_type in [:full, :semi]
            out[:forecastobs]    = DSGE.transplant_forecast_observables(histstates, forecastobs, system, T)
            out[:forecastpseudo] = DSGE.transplant_forecast(histpseudo, forecastpseudo, T)
        else
            out[:forecastobs]    = forecastobs
            out[:forecastpseudo] = forecastpseudo
        end

        # Keep only forecasted periods between keep_startdate and keep_enddate
        forecast_dates = DSGE.quarter_range(date_forecast_start(m), date_forecast_end(m))
        keep_inds = keep_startdate .<= forecast_dates .<= keep_enddate
        out[:forecastobs]    = out[:forecastobs][:, keep_inds]
        out[:forecastpseudo] = out[:forecastpseudo][:, keep_inds]
    end

    # Compute trend, dettrend, and shockdecs
    if :shockdec in outputs
        m <= Setting(:shockdec_startdate, Nullable(keep_startdate))
        m <= Setting(:shockdec_enddate,   Nullable(keep_enddate))

        _, out[:trendobs],    out[:trendpseudo]    = trends(system)
        _, out[:dettrendobs], out[:dettrendpseudo] = deterministic_trends(m, system, s_0)

        # Applying ϵ_{1:T-k}
        t0 = date_mainsample_start(m)
        t1 = shockdec_splitdate
        _, out[:revisionobs], out[:revisionpseudo] =
            DSGE.shock_decompositions(m, system, histshocks, shock_start_date = t0, shock_end_date = t1)

        # Applying ϵ_{T-k+1:end}
        t0 = DSGE.iterate_quarters(shockdec_splitdate, 1)
        t1 = cond_type == :none ? date_mainsample_end(m) : date_conditional_end(m)
        _, out[:newsobs], out[:newspseudo] =
            DSGE.shock_decompositions(m, system, histshocks, shock_start_date = t0, shock_end_date = t1)

        # Shockdec output starts as size Ny x Nh x Ne
        # Sum over shock dimension to get size Ny x Nh
        for output_var in [:revisionobs, :revisionpseudo, :newsobs, :newspseudo]
            out[output_var] = squeeze(sum(out[output_var], 3), 3)
        end
    end

    # Return
    return out
end