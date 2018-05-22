function class_measurement_matrices(sys::System, class::Symbol)
    ZZ, DD = if class == :obs
        sys[:ZZ], sys[:DD]
    elseif class == :pseudo
        sys[:ZZ_pseudo], sys[:DD_pseudo]
    elseif class == :states
        Ns = size(sys[:TTT], 1)
        eye(Ns), zeros(Ns)
    else
        error("Invalid class: $class. Must be :obs, :pseudo, or :state")
    end
    return ZZ, DD
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
            return false
        end
    end
    return true
end

function check_states_shocks_decomp(sys_new::System, s_tgt::Matrix{Float64}, s_tgT::Matrix{Float64},
                                    classes::Vector{Symbol}, k::Int, h::Int,
                                    state_comps::Dict{Symbol, Vector{Float64}},
                                    shock_comps::Dict{Symbol, Vector{Float64}};
                                    atol::Float64 = 1e-8)
    TTT, CCC = sys_new[:TTT], sys_new[:CCC]

    # s_{T-k+h|T}
    s_Tmkph_T = if h <= k
        s_tgT[:, end-k+h] + CCC # history
    else
        TTT^(h-k) * s_tgt[:, end] + CCC # forecast
    end

    # s_{T-k+h|T-k}
    s_Tmkph_Tmk = TTT^h * s_tgt[:, end-k] + CCC

    for class in classes
        ZZ, DD = class_measurement_matrices(sys_new, class)
        y_Tmkph_T   = ZZ*s_Tmkph_T   + DD
        y_Tmkph_Tmk = ZZ*s_Tmkph_Tmk + DD

        # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
        state_shock_comp = y_Tmkph_T - y_Tmkph_Tmk

        if !isapprox(state_shock_comp, state_comps[class] + shock_comps[class], atol = atol)
            @show state_comps[class] + shock_comps[class]
            @show state_shock_comp
            return false
        end
    end
    return true
end