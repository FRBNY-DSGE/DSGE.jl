function compute_system(m; use_expected_rate_data=false)

# Solve model
TTT, RRR, CCC = solve(m)

# Index out non-anticipated shocks rows and columns of system matrices if desired
if !use_expected_rate_data
    n_states     = DSGE.n_states(m)
    n_ant        = n_anticipated_shocks(m)
    n_states_aug = DSGE.n_states_augmented(m)
    n_exo        = n_shocks_exogenous(m)

    state_inds = [1:(n_states-n_ant); (n_states+1):n_states_aug]
    shock_inds = 1:(n_exo-n_ant)

    TTT = TTT[state_inds, state_inds]
    RRR = RRR[state_inds, shock_inds]
    CCC = CCC[state_inds, :]
end

transition_equation = Transition(TTT, RRR, CCC)

# Solve measurement equation
measurement_equation = measurement(m, TTT, RRR, CCC; shocks = use_expected_rate_data)

return System(transition_equation, measurement_equation)

end
