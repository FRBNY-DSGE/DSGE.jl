function shock_loading(m::RealBond, TTT_jump::Matrix{Float64})
    exo = m.exogenous_shocks
    endo = m.endogenous_states

    _RRR = zeros(n_backward_looking_states(m), n_shocks_exogenous(m))
    RRR  = zeros(n_model_states(m), n_shocks_exogenous(m))

    # Technology shock loading on states
    _RRR[endo[:z′_t], exo[:z_sh]] .= 1.

    # Monetary policy shock loading on states
    _RRR[endo[:mon′_t], exo[:mon_sh]] .= 1.

    # Loading on states and jumps
    RRR[1:n_backward_looking_states(m), :] = _RRR
    RRR[n_backward_looking_states(m)+1:end, :] = TTT_jump*_RRR

    return RRR
end
