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