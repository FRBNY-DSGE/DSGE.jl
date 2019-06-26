function simulate_irfs(T::Matrix{Float64}, R::Matrix{Float64}, vars_SS::OrderedDict{Symbol, Any},
                       periods::Int, steps::Int, transformation::Matrix{Float64};
                       DisplayLev::Int = 2)

    DisplayLev >= 1 && println("Simulating IRFs...")
    # Set up shock
    dt = periods/steps
    agg_shock = zeros(1, steps)
    agg_shock[1] = 1/sqrt(dt)

    # Simulate
    simulated_states = simulate(T, R, periods, steps, agg_shock; method = :implicit)
    simulated_states = transformation * simulated_states
    # Get sizes
    I, J = size(vars_SS[:V_SS])
    n_v = I*J + 1
    n_g = I*J

    # Normalize if necessary
    inflation      = simulated_states[n_v,       :]
    monetary_shock = simulated_states[n_v+n_g,   :]
    wage           = simulated_states[n_v+n_g+1, :] / vars_SS[:w_SS]
    consumption    = simulated_states[n_v+n_g+2, :] / vars_SS[:C_SS]
    labor_supply   = simulated_states[n_v+n_g+3, :] / vars_SS[:N_SS]
    output         = simulated_states[n_v+n_g+4, :] / vars_SS[:Y_SS]

    return inflation, monetary_shock, wage, consumption, labor_supply, output
end