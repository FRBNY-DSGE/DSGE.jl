function decompose_forecast_new(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                                params_new::Vector{Float64}, params_old::Vector{Float64},
                                cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol},
                                hs::UnitRange{Int};
                                individual_shocks::Bool = false,
                                check::Bool = false, atol::Float64 = 1e-8) where M<:AbstractModel
    # Compute numbers of periods
    T0, T, k, T1_new, T1_old, _ = DSGE.decomposition_periods(m_new, m_old, cond_new, cond_old)
    @assert size(df_new, 1) == T0 + T + T1_new
    @assert size(df_old, 1) == T0 + T - k + T1_old

    # Forecast
    for m in [m_new, m_old]
        m <= Setting(:shockdec_startdate, Nullable(date_forecast_start(m_new)))
        m <= Setting(:shockdec_enddate,   Nullable(date_forecast_end(m_old)))
    end
    output_vars = vcat([:shockdecobs, :dettrendobs, :trendobs], check ? [:forecastobs] : [])
    out1 = DSGE.forecast_one_draw(m_new, :mode, cond_new, output_vars, params_new, df_new) # new data, new params
    out2 = DSGE.forecast_one_draw(m_old, :mode, cond_old, output_vars, params_new, df_old) # old data, new params
    out3 = DSGE.forecast_one_draw(m_old, :mode, cond_old, output_vars, params_old, df_old) # old data, old params

    # out[:shockdecobs] starts as size Ny x Nh x Ne
    # Sum over shock dimension to get size Ny x Nh
    for out in [out1, out2, out3]
        out[:shockdecobs] = squeeze(sum(out[:shockdecobs], 3), 3)
    end

    # Decomposition
    # out[:dettrendobs], out[:shockdecobs], and out[:trendobs] have size Ny x Nh
    # out[:dettrendobs][:, h] is Z T^{T+h} s_{0|T} + D
    # out[:shockdecobs][:, h] is Z \sum_{t=1}^T T^{T+h-t} R ϵ_{t|T} + D
    # out[:trendobs][:, h]    is Z*D

    # Data revision + news
    data_news_comp = (out1[:dettrendobs] - out2[:dettrendobs]) + (out1[:shockdecobs] - out2[:shockdecobs])

    # Parameter re-estimation
    y_θ_new = out2[:dettrendobs] + out2[:shockdecobs] .+ out2[:trendobs]
    y_θ_old = out3[:dettrendobs] + out3[:shockdecobs] .+ out3[:trendobs]
    para_comp = y_θ_new - y_θ_old

    # Check decomposition
    if check
        Nh = size(out1[:dettrendobs], 2)
        y_new = out1[:forecastobs][:, 1:Nh]
        y_old = out3[:forecastobs][:, k+1:k+Nh]
        @assert para_comp + data_news_comp ≈ y_new - y_old
    end

    return data_news_comp, para_comp
end