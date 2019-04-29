function shock_loading(m::HetDSGE, TTT_jump::Matrix{Float64})
    exo = m.exogenous_shocks
    endo = m.endogenous_states

    _RRR = zeros(n_backward_looking_states(m), n_shocks_exogenous(m))
    RRR  = zeros(n_model_states(m), n_shocks_exogenous(m))

    _RRR[endo[:b′_t], exo[:b_sh]] .= 1.
    _RRR[endo[:g′_t], exo[:g_sh]] .= 1.
    _RRR[endo[:z′_t], exo[:z_sh]] .= 1.
    _RRR[endo[:μ′_t], exo[:μ_sh]] .= 1.
    _RRR[endo[:λ_w′_t], exo[:λ_w_sh]] .= 1.
    _RRR[endo[:λ_f′_t], exo[:λ_f_sh]] .= 1.
    _RRR[endo[:rm′_t], exo[:rm_sh]] .= 1.

    # Loading on states and jumps
    RRR[1:n_backward_looking_states(m), :] = _RRR
    RRR[n_backward_looking_states(m)+1:end, :] = TTT_jump*_RRR

    return RRR
end
