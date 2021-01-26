function augment_states(m::SmetsWouters{T}, TTT::Matrix{T}, RRR::Matrix{T},
                        CCC::Vector{T};
                        regime_switching::Bool = false,
                        reg::Int = 1) where {T<:AbstractFloat}
    endo = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    exo = m.exogenous_shocks

    n_endo = n_states(m)
    n_exo = n_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert n_endo == length(CCC)

    # Initialize augmented matrices
    n_addl_states = length(endo_addl)
    n_addl_eqs = n_addl_states
    TTT_aug = zeros(n_endo + n_addl_eqs, n_endo + n_addl_states)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(n_addl_eqs, n_exo)]
    CCC_aug = [CCC; zeros(n_addl_eqs)]

    for p in m.parameters
        if !isempty(p.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end


    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl[:y_t1], endo[:y_t]] = 1.0
    TTT_aug[endo_addl[:c_t1], endo[:c_t]] = 1.0
    TTT_aug[endo_addl[:i_t1], endo[:i_t]] = 1.0
    TTT_aug[endo_addl[:w_t1], endo[:w_t]] = 1.0
    TTT_aug[endo_addl[:π_t1], endo[:π_t]] = 1.0
    TTT_aug[endo_addl[:L_t1], endo[:L_t]]  = 1.0

    ## We construct state for expected inflation using the fact
    ##
    ##   E_t[s_{t+1}] = CCC + TTT*s_{t} = CCC + TTT * (CCC + TTT*s_{t-1} + RRR)
    ##                = (CCC+TTT*CCC) + (TTT^2)s_{t-1} + TTT*RRR
    ##
    ## So to construct the state for E_t[p_{t+1}], we need to make all these
    ## objects, and just index out the row relevant to pi_t

    T2 = TTT^2
    TR = TTT*RRR
    CTC = CCC+TTT*CCC


    TTT_aug[endo_addl[:Et_π_t],:] = [T2[endo[:π_t],:]; zeros(n_addl_states)]

    RRR_aug[endo_addl[:Et_π_t],:] = TR[endo[:π_t],:]
    CCC_aug[endo_addl[:Et_π_t],:] = CTC[endo[:π_t],:]

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return TTT_aug, RRR_aug, CCC_aug
end
