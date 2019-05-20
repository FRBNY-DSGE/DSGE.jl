function augment_states(m::RepDSGEGovDebt{T}, TTT::Matrix{T},
                        RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}

    endo     = m.endogenous_states # augment_model_states(m.endogenous_states, n_model_states(m))
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks

    n_endo = n_model_states(m)
    n_exo  = n_shocks_exogenous(m)

    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert n_endo == length(CCC)

    # Initialize augmented matrices
    n_new_states = length(endo_new)
    n_new_eqs = n_new_states
    TTT_aug = zeros(n_endo + n_new_eqs, n_endo + n_new_states)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(n_new_eqs, n_exo)]
    CCC_aug = [CCC; zeros(n_new_eqs)]

    ### TTT modifications

    # Track Lags
    TTT_aug[endo_new[:c_t1], endo[:l′_t]] = -1.0

    return TTT_aug, RRR_aug, CCC_aug
end
