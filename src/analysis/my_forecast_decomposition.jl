# TODO:
# - Handle conditioning
# - Compute multiple h at once
# - Handle pseudo-observables
# - Maybe: break out checking into own functions
# - Maybe: remove duplicate filter runs

function decompose_forecast(m_new::AbstractModel, m_old::AbstractModel,
                            params_new::Vector{Float64}, params_old::Vector{Float64},
                            df_new::DataFrame, df_old::DataFrame,
                            k::Int, h::Int;
                            check::Bool = false)
    # Update parameters
    prepare_decomposition!(m_new, m_old, params_new, params_old, k)

    # Decompose into components
    state_comp, shock_comp = decompose_states_shocks(m_new, df_new, k, h, check = check)
    data_comp = decompose_data_revisions(m_new, df_new, df_old, k, h, check = check)
    param_comp = decompose_param_reest(m_new, m_old, df_old, k, h, check = check)

    # Compute total difference
    total = state_comp + shock_comp + data_comp + param_comp

    if check
        T  = n_mainsample_periods(m_new)

        # y^{old,old}_{T-k+h|T-k}
        sys_old = compute_system(m_old)
        kal_old = DSGE.filter(m_old, df_old, sys_old)
        _, y_old, _ = forecast(m_old, sys_old, kal_old[:zend]) # y^{old,old}_{T-k+j|T-k}, j = 1:H
        y_old_Tmkph_Tmk = y_old[:, h]

        # y^{new,new}_{T-k+h|T}
        sys_new = compute_system(m_new)
        kal_new = DSGE.filter(m_new, df_new, sys_new)
        y_new_Tmkph_T = if h <= k
            # New history vs. old forecast
            s_new, _ = smooth(m_new, df_new, sys_new, kal_new, draw_states = false)
            sys_new[:ZZ] * s_new[:, T-k+h] + sys_new[:DD] # y^{new,new}_{t|T}, t = 1:T
        else
            # New forecast vs. old forecast
            _, y_new, _ = forecast(m_new, sys_new, kal_new[:zend]) # y^{new,new}_{T+j|T}, j = 1:H
            y_new[:, h-k]
        end

        # Total difference = y^{new,new}_{T+h-k|T} - y^{old,old}_{T-k+h|T-k}
        exp_total = y_new_Tmkph_T - y_old_Tmkph_Tmk
        try
            @assert total ≈ exp_total
        catch ex
            @show total
            @show exp_total
            throw(ex)
        end
    end

    # Put into DataFrame and return
    out = DataFrame(var = collect(keys(m_new.observables)), total = total,
                    state = state_comp, shock = shock_comp, data = data_comp, param = param_comp)
    return out
end

function prepare_decomposition!(m_new::AbstractModel, m_old::AbstractModel,
                                params_new::Vector{Float64}, params_old::Vector{Float64},
                                k::Int)
    # Check models well-formed
    @assert typeof(m_new) == typeof(m_old)
    @assert DSGE.subtract_quarters(date_forecast_start(m_new), date_forecast_start(m_old)) == k

    # Update parameters
    DSGE.update!(m_new, params_new)
    DSGE.update!(m_old, params_old)
end

function decompose_states_shocks(m_new::AbstractModel, df_new::DataFrame, k::Int, h::Int;
                                 check::Bool = false)
    # New parameters, new data
    # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
    #     = Z [ T^h (s_{T-k|T} - s_{T-k|T-k}) + sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} ) ]
    #     = state component + shock component
    T0 = n_presample_periods(m_new)
    T  = n_mainsample_periods(m_new)
    sys_new = compute_system(m_new)
    ZZ, DD, TTT, RRR, CCC = sys_new[:ZZ], sys_new[:DD], sys_new[:TTT], sys_new[:RRR], sys_new[:CCC]

    kal_new      = DSGE.filter(m_new, df_new, sys_new)
    s_tgt        = kal_new[:filt][:, (T0+1):end] # s_{t|t}
    s_tgT, ϵ_tgT = smooth(m_new, df_new, sys_new, kal_new, draw_states = false) # s_{t|T}, ϵ_{t|T}

    # State component = Z T^h (s_{T-k|T} - s_{T-k|T-k})
    s_Tmk_T   = s_tgT[:, T-k] # s_{T-k|T}
    s_Tmk_Tmk = s_tgt[:, T-k] # s_{T-k|T-k}
    state_comp = ZZ * TTT^h * (s_Tmk_T - s_Tmk_Tmk)

    # Shock component = Z sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} )
    summation_terms = map(j -> TTT^(h-j) * RRR * ϵ_tgT[:, T-k+j], 1:min(k, h))
    shock_comp = ZZ * sum(summation_terms)

    # Return
    if check
        @assert size(s_tgt, 2) == T
        @assert size(s_tgT, 2) == T
        @assert size(ϵ_tgT, 2) == T
        @assert s_tgt[:, T] ≈ s_tgT[:, T]

        y_Tmkph_T = if h <= k
            ZZ * (s_tgT[:, T-k+h] + CCC) + DD # smoothed history
        else
            ZZ * (TTT^(h-k) * s_tgt[:, T] + CCC) + DD # forecast
        end
        y_Tmkph_Tmk = ZZ * (TTT^h * s_tgt[:, T-k] + CCC) + DD
        state_shock_comp = y_Tmkph_T - y_Tmkph_Tmk

        try
            @assert state_comp + shock_comp ≈ state_shock_comp
        catch ex
            @show state_comp + shock_comp
            @show state_shock_comp
        end
    end
    return state_comp, shock_comp
end

function decompose_data_revisions(m_new::AbstractModel, df_new::DataFrame, df_old::DataFrame,
                                  k::Int, h::Int; check::Bool = false)
    # New parameters, new and old data
    # Data revision component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}
    #     = Z T^h ( s^{new}_{T-k|T-k} - s^{old}_{T-k|T-k} )
    T = n_mainsample_periods(m_new)
    sys_new = compute_system(m_new)
    ZZ, TTT = sys_new[:ZZ], sys_new[:TTT]

    kal_new = DSGE.filter(m_new, df_new, sys_new, include_presample = false)
    s_new_tgt = kal_new[:filt] # s^{new}_{t|t}
    check && @assert size(s_new_tgt, 2) == T

    kal_old = DSGE.filter(m_new, df_old, sys_new, include_presample = false)
    s_old_tgt = kal_old[:filt] # s^{old}_{t|t}
    check && @assert size(s_old_tgt, 2) == T-k

    # Return
    data_comp = ZZ * TTT^h * (s_new_tgt[:, T-k] - s_old_tgt[:, T-k])
    return data_comp
end

function decompose_param_reest(m_new::AbstractModel, m_old::AbstractModel, df_old::DataFrame,
                               k::Int, h::Int; check::Bool = false)
    # New and old parameters, old data
    # Parameter re-estimation component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}

    # y^{new}_{T-k+h|T-k} = Z^{new} (T^{new}^h s^{new}_{T-k|T-k} + C^{new}) + D^{new}
    T = n_mainsample_periods(m_new)
    sys_new = compute_system(m_new)
    ZZ_new, DD_new, TTT_new, CCC_new = sys_new[:ZZ], sys_new[:DD], sys_new[:TTT], sys_new[:CCC]

    kal_new = DSGE.filter(m_new, df_old, sys_new, include_presample = false) # s^{new}_{t|t}
    s_new_tgt = kal_new[:filt]
    check && @assert size(s_new_tgt, 2) == T-k
    s_new_Tmk_Tmk = s_new_tgt[:, T-k]
    y_new_Tmkph_Tmk = ZZ_new * (TTT_new^h * s_new_Tmk_Tmk + CCC_new) + DD_new

    # y^{old}_{T-k+h|T-k} = Z^{old} (T^{old}^h s^{old}_{T-k|T-k} + C^{old}) + D^{old}
    sys_old = compute_system(m_old)
    ZZ_old, DD_old, TTT_old, CCC_old = sys_old[:ZZ], sys_old[:DD], sys_old[:TTT], sys_old[:CCC]

    kal_old = DSGE.filter(m_old, df_old, sys_old, include_presample = false) # s^{old}_{t|t}
    s_old_tgt = kal_old[:filt]
    check && @assert size(s_old_tgt, 2) == T-k
    s_old_Tmk_Tmk = s_old_tgt[:, T-k]
    y_old_Tmkph_Tmk = ZZ_old * (TTT_old^h * s_old_Tmk_Tmk + CCC_old) + DD_old

    # Return
    param_comp = y_new_Tmkph_Tmk - y_old_Tmkph_Tmk
    return param_comp
end