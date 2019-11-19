# Set up model in canonical form

function f{T<:Real}(x::Vector{T}; vars_ss = deepcopy(vars_SS))
    v_residual = equilibrium_conditions(x, vars_ss, grids, params, n_v, n_g, n_p, niter_hours = 10)
    return vcat(values(v_residual)...)
end
derivs = ForwardDiff.jacobian(f, x)
Γ1 = -derivs[:,1:n_vars];
Γ0 = derivs[:,n_vars+1:2*n_vars];
Π = -derivs[:,2*n_vars+1:2*n_vars+n_exp_errors];
Ψ = -derivs[:,2*n_vars+n_exp_errors+1:2*n_vars+n_exp_errors+n_shocks];
C = zeros(n_vars, 1)
