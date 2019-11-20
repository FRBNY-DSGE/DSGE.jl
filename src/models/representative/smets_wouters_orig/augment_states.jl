function augment_states(m::SmetsWoutersOrig{T}, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where T<:AbstractFloat
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

    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl[:y_t1], endo[:y_t]] = 1.
    TTT_aug[endo_addl[:c_t1], endo[:c_t]] = 1.
    TTT_aug[endo_addl[:i_t1], endo[:i_t]] = 1.
    TTT_aug[endo_addl[:w_t1], endo[:w_t]] = 1.

    return TTT_aug, RRR_aug, CCC_aug
end
