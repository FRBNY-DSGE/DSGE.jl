function class_system_matrices(sys::System, class::Symbol)
    TTT, RRR, CCC = sys[:TTT], sys[:RRR], sys[:CCC]
    ZZ, DD = if class == :obs
        sys[:ZZ], sys[:DD]
    elseif class == :pseudo
        sys[:ZZ_pseudo], sys[:DD_pseudo]
    elseif class == :states
        I, zeros(size(TTT, 1))
    else
        error("Invalid class: $class. Must be :obs, :pseudo, or :state")
    end
    return ZZ, DD, TTT, RRR, CCC
end

function compute_history_and_forecast(m::AbstractModel, df::DataFrame, class::Symbol;
                                      cond_type::Symbol = :none)
    sys = compute_system(m)
    ZZ, DD, _, _, _ = class_system_matrices(sys, class)

    histstates, _ = smooth(m, df, sys, cond_type = cond_type, draw_states = false)
    hist = ZZ * histstates .+ DD

    fcast = Dict{Symbol, Matrix{Float64}}()
    fcast[:states], fcast[:obs], fcast[:pseudo], _ =
        forecast(m, sys, histstates[:, end], cond_type = cond_type)

    return hcat(hist, fcast[class])
end

function check_total_decomp(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                            cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol},
                            decomp::Dict{Symbol, Array{Float64}};
                            hs = 1:forecast_horizons(m_old),
                            atol::Float64 = 1e-8) where M<:AbstractModel
    correct = true
    for class in classes
        # y^{new,new}_{t|T}, t = 1:T+H
        y_new = compute_history_and_forecast(m_new, df_new, class, cond_type = cond_new)

        # y^{old,old}_{t|T-k}, t = 1:T-k+H
        y_old = compute_history_and_forecast(m_old, df_old, class, cond_type = cond_old)

        # Total difference = y^{new,new}_{T-k+h|T} - y^{old,old}_{T-k+h|T-k}
        y_new_Tmkph_T   = y_new[:, T-k+hs]
        y_old_Tmkph_Tmk = y_old[:, T-k+hs]
        exp_total = y_new_Tmkph_T - y_old_Tmkph_Tmk

        # Compare sum of decomposed differences to expected total difference
        if !isapprox(exp_total, decomp[Symbol(:decomptotal, class)], atol = atol)
            @show decomp[Symbol(:decomptotal, class)]
            @show exp_total
            correct = false
        end
    end
    return correct
end

function check_states_shocks_decomp(sys_new::System, s_tgt::Matrix{Float64}, s_tgT::Matrix{Float64},
                                    class::Symbol, k::Int, h::Int,
                                    state_comp::Vector{Float64}, shock_comp::Vector{Float64};
                                    atol::Float64 = 1e-8)
    ZZ, DD, TTT, RRR, CCC = class_system_matrices(sys_new, class)

    # y_{T-k+h|T}
    y_Tmkph_T = if h <= k
        ZZ * (s_tgT[:, end-k+h] + CCC) + DD # history
    else
        ZZ * (TTT^(h-k) * s_tgt[:, end] + CCC) + DD # forecast
    end

    # y_{T-k+h|T-k}
    y_Tmkph_Tmk = ZZ * (TTT^h * s_tgt[:, end-k] + CCC) + DD

    # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
    state_shock_comp = y_Tmkph_T - y_Tmkph_Tmk

    if isapprox(state_shock_comp, state_comp + shock_comp, atol = atol)
        return true
    else
        @show state_comp + shock_comp
        @show state_shock_comp
        return false
    end
end