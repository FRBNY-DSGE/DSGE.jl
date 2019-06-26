
##This is the model file for the Aiyagari model

params, approx_params, grids, approx_valuef, I, J = set_parameters()
println("Parameters set.")

# n_v = Int64(I*J) # The number of jump variables (value function)
# n_g = Int64(I*J)     # The number of endogenous state variables (distribution)
# n_p = 2              # The number of static relations: bond-market clearing, labor
                     # market clearing, consumption, output, total assets
# n_exp_errors = n_v
# n_vars = n_v + n_g + n_p

# compute steady state
println("computing steady state...")
vars_SS = compute_steady_state(grids, params, approx_params; DisplayLev = 0)

