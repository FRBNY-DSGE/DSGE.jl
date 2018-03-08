# TODO:
# - Handle pseudo-observables
# - Maybe: break out checking into own functions
# - Save decomposition output
# - Plotting

function decompose_forecast(m_new::AbstractModel, m_old::AbstractModel,
                            df_new::DataFrame, df_old::DataFrame,
                            params_new::Vector{Float64}, params_old::Vector{Float64};
                            hs = 1:forecast_horizons(m_old),
                            cond_new::Symbol = :none, cond_old::Symbol = cond_new,
                            check::Bool = false, atol::Float64 = 1e-8)

    # Compute numbers of periods
    # h_cond and k_cond are h and k adjusted for (differences in) conditioning
    T0, T, k, T1_new, T1_old, k_cond =
        decomposition_periods(m_new, m_old, cond_new = cond_new, cond_old = cond_old)

    if check
        @assert size(df_new, 1) == T0 + T + T1_new
        @assert size(df_old, 1) == T0 + T + T1_old - k
    end

    # Update parameters, filter, and smooth
    sys_new, sys_old, kal_new_new, kal_new_old, kal_old_old, s_tgT, ϵ_tgT =
        prepare_decomposition!(m_new, m_old, df_new, df_old, params_new, params_old,
                               cond_new = cond_new, cond_old = cond_old)

    # Initialize DataFrames
    dfs = DataStructures.OrderedDict{Symbol, DataFrame}()
    for var in keys(m_new.observables)
        dfs[var] = DataFrame(h = Int[], date = Date[], total = Float64[],
                             state = Float64[], shock = Float64[], data = Float64[], param = Float64[])
    end

    # Compute whole sequence of historical and forecasted observables if desired
    if check
        # y^{new,new}_{t|T-k}, t = 1:T:H
        y_new = compute_history_and_forecast(m_new, df_new, cond_type = cond_new)

        # y^{old,old}_{t|T-k}, t = 1:T-k+H
        y_old = compute_history_and_forecast(m_old, df_old, cond_type = cond_old)
    end

    for h in hs
        # Adjust h for (differences in) conditioning
        h_cond = h - T1_old
        @assert h_cond >= 0

        # Decompose into components
        state_comp, shock_comp = decompose_states_shocks(sys_new, kal_new_new, s_tgT, ϵ_tgT,
                                                         T0, k_cond, h_cond, check = check)
        data_comp = decompose_data_revisions(sys_new, kal_new_new, kal_new_old,
                                             T0, k_cond, h_cond, check = check)
        param_comp = decompose_param_reest(sys_new, sys_old, kal_new_old, kal_old_old,
                                           T0, k_cond, h_cond, check = check)

        # Compute total difference
        total = state_comp + shock_comp + data_comp + param_comp

        # Add to DataFrames
        date = DSGE.iterate_quarters(date_mainsample_end(m_old), h)
        for (var, i) in m_new.observables
            push!(dfs[var], [h, date, total[i], state_comp[i], shock_comp[i], data_comp[i], param_comp[i]])
        end

        if check
            # Total difference = y^{new,new}_{T+h-k|T} - y^{old,old}_{T-k+h|T-k}
            y_new_Tmkph_T   = y_new[:, T-k+h]
            y_old_Tmkph_Tmk = y_old[:, T-k+h]
            exp_total = y_new_Tmkph_T - y_old_Tmkph_Tmk
            try
                @assert isapprox(exp_total, total, atol = atol)
            catch ex
                @show total
                @show exp_total
                throw(ex)
            end
        end
    end

    return dfs
end

function decomposition_periods(m_new::AbstractModel, m_old::AbstractModel;
                               cond_new::Symbol = :none, cond_old::Symbol = cond_old)
    # Number of presample periods T0 must be the same
    T0 = n_presample_periods(m_new)
    @assert n_presample_periods(m_old) == T0

    # New model has T main-sample periods
    # Old model has T-k main-sample periods
    T = n_mainsample_periods(m_new)
    k = DSGE.subtract_quarters(date_forecast_start(m_new), date_forecast_start(m_old))
    @assert k >= 0

    # Number of conditional periods T1 may differ
    T1_new = cond_new == :none ? 0 : n_conditional_periods(m_new)
    T1_old = cond_old == :none ? 0 : n_conditional_periods(m_old)

    # Adjust k for (differences in) conditioning
    k_cond = k + T1_new - T1_old

    return T0, T, k, T1_new, T1_old, k_cond
end

function prepare_decomposition!(m_new::AbstractModel, m_old::AbstractModel,
                                df_new::DataFrame, df_old::DataFrame,
                                params_new::Vector{Float64}, params_old::Vector{Float64};
                                cond_new::Symbol = :none, cond_old::Symbol = cond_old)
    # Check models well-formed
    @assert typeof(m_new) == typeof(m_old)

    # Update parameters
    DSGE.update!(m_new, params_new)
    DSGE.update!(m_old, params_old)
    sys_new = compute_system(m_new)
    sys_old = compute_system(m_old)

    # Filter and smooth
    kal_new_new = DSGE.filter(m_new, df_new, sys_new, cond_type = cond_new)
    kal_new_old = DSGE.filter(m_new, df_old, sys_new, cond_type = cond_old)
    kal_old_old = DSGE.filter(m_old, df_old, sys_old, cond_type = cond_old)
    s_tgT, ϵ_tgT = smooth(m_new, df_new, sys_new, kal_new_new,
                          cond_type = cond_new, draw_states = false)

    return sys_new, sys_old, kal_new_new, kal_new_old, kal_old_old, s_tgT, ϵ_tgT
end

function decompose_states_shocks(sys_new::System, kal_new::DSGE.Kalman,
                                 s_tgT::Matrix{Float64}, ϵ_tgT::Matrix{Float64},
                                 T0::Int, k::Int, h::Int;
                                 check::Bool = false, atol::Float64 = 1e-8)
    # New parameters, new data
    # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
    #     = Z [ T^h (s_{T-k|T} - s_{T-k|T-k}) + sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} ) ]
    #     = state component + shock component
    ZZ, DD, TTT, RRR, CCC = sys_new[:ZZ], sys_new[:DD], sys_new[:TTT], sys_new[:RRR], sys_new[:CCC]
    s_tgt = kal_new[:filt][:, (T0+1):end] # s_{t|t}

    # State component = Z T^h (s_{T-k|T} - s_{T-k|T-k})
    s_Tmk_T   = s_tgT[:, end-k] # s_{T-k|T}
    s_Tmk_Tmk = s_tgt[:, end-k] # s_{T-k|T-k}
    state_comp = ZZ * TTT^h * (s_Tmk_T - s_Tmk_Tmk)

    # Shock component = Z sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} )
    shock_sum = zeros(size(s_Tmk_T))
    for j = 1:min(k, h)
        shock_sum += TTT^(h-j) * RRR * ϵ_tgT[:, end-k+j]
    end
    shock_comp = ZZ * shock_sum

    # Return
    if check
        @assert size(s_tgt, 2) == size(s_tgT, 2) == size(ϵ_tgT, 2)
        @assert isapprox(s_tgt[:, end], s_tgT[:, end], atol = atol)

        y_Tmkph_T = if h <= k
            ZZ * (s_tgT[:, end-k+h] + CCC) + DD # history
        else
            ZZ * (TTT^(h-k) * s_tgt[:, end] + CCC) + DD # forecast
        end
        y_Tmkph_Tmk = ZZ * (TTT^h * s_tgt[:, end-k] + CCC) + DD
        state_shock_comp = y_Tmkph_T - y_Tmkph_Tmk

        try
            @assert isapprox(state_shock_comp, state_comp + shock_comp, atol = atol)
        catch ex
            @show state_comp + shock_comp
            @show state_shock_comp
            throw(ex)
        end
    end
    return state_comp, shock_comp
end

function decompose_data_revisions(sys_new::System, kal_new::DSGE.Kalman, kal_old::DSGE.Kalman,
                                  T0::Int, k::Int, h::Int;
                                  check::Bool = false, atol::Float64 = 1e-8)
    # New parameters, new and old data
    # Data revision component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}
    #     = Z T^h ( s^{new}_{T-k|T-k} - s^{old}_{T-k|T-k} )
    ZZ, TTT = sys_new[:ZZ], sys_new[:TTT]

    s_new_tgt = kal_new[:filt][:, (T0+1):end] # s^{new}_{t|t}
    s_old_tgt = kal_old[:filt][:, (T0+1):end] # s^{old}_{t|t}
    try
        check && @assert size(s_new_tgt, 2) == size(s_old_tgt, 2) + k
    catch ex
        @show size(s_new_tgt, 2)
        @show size(s_old_tgt, 2)
        @show k
        throw(ex)
    end

    # Return
    data_comp = ZZ * TTT^h * (s_new_tgt[:, end-k] - s_old_tgt[:, end-k])
    return data_comp
end

function decompose_param_reest(sys_new::System, sys_old::System,
                               kal_new::DSGE.Kalman, kal_old::DSGE.Kalman,
                               T0::Int, k::Int, h::Int;
                               check::Bool = false, atol::Float64 = 1e-8)
    # New and old parameters, old data
    # Parameter re-estimation component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}
    ZZ_new, DD_new, TTT_new, CCC_new = sys_new[:ZZ], sys_new[:DD], sys_new[:TTT], sys_new[:CCC]
    ZZ_old, DD_old, TTT_old, CCC_old = sys_old[:ZZ], sys_old[:DD], sys_old[:TTT], sys_old[:CCC]

    s_new_tgt = kal_new[:filt][:, (T0+1):end] # s^{new}_{t|t}
    s_old_tgt = kal_old[:filt][:, (T0+1):end] # s^{old}_{t|t}
    check && @assert size(s_new_tgt, 2) == size(s_old_tgt, 2)

    # y^{new}_{T-k+h|T-k} = Z^{new} (T^{new}^h s^{new}_{T-k|T-k} + C^{new}) + D^{new}
    s_new_Tmk_Tmk = s_new_tgt[:, end-k]
    y_new_Tmkph_Tmk = ZZ_new * (TTT_new^h * s_new_Tmk_Tmk + CCC_new) + DD_new

    # y^{old}_{T-k+h|T-k} = Z^{old} (T^{old}^h s^{old}_{T-k|T-k} + C^{old}) + D^{old}
    s_old_Tmk_Tmk = s_old_tgt[:, end]
    y_old_Tmkph_Tmk = ZZ_old * (TTT_old^h * s_old_Tmk_Tmk + CCC_old) + DD_old

    # Return
    param_comp = y_new_Tmkph_Tmk - y_old_Tmkph_Tmk
    return param_comp
end



### HELPER FUNCTIONS

function Base.getindex(kal::DSGE.Kalman, inds::Union{Int, UnitRange{Int}})
    t0 = first(inds)
    t1 = last(inds)

    return DSGE.Kalman(sum(kal[:marginal_L][inds]),  # L
                       kal[:filt][:, t1],            # zend
                       kal[:vfilt][:, :, t1],        # Pend
                       kal[:pred][:, inds],          # pred
                       kal[:vpred][:, :, inds],      # vpred
                       kal[:yprederror][:, inds],    # yprederror
                       kal[:ystdprederror][:, inds], # ystdprederror
                       sqrt.(mean((kal[:yprederror][:, inds].^2)', 1)), # rmse
                       sqrt.(mean((kal[:ystdprederror][:, inds].^2)', 1)), # rmsd
                       kal[:filt][:, inds],          # filt
                       kal[:vfilt][:, :, inds],      # vfilt
                       kal[:filt][:, t0],            # z0
                       kal[:vfilt][:, :, t0],        # P0
                       kal[:marginal_L][inds])       # marginal_L
end

function compute_history_and_forecast(m::AbstractModel, df::DataFrame; cond_type::Symbol = :none)
    sys = compute_system(m)
    kal = DSGE.filter(m, df, sys, cond_type = cond_type)
    histstates, _ = smooth(m, df, sys, kal, cond_type = cond_type, draw_states = false)
    hist = sys[:ZZ] * histstates .+ sys[:DD]
    _, fcast, _ = forecast(m, sys, kal[:zend])
    return hcat(hist, fcast)
end
