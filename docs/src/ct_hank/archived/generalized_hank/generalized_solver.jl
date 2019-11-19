
#We need to let the user speciy periods and steps for the irfs

function generalized_solver(canonical_form::Dict{Symbol, Any}, params::Dict{Symbol, Float64},
                            init_params::Dict{Symbol, Float64}, grids::Dict{Symbol, Any},
                            approx_params::Dict{Symbol, Any}, red_params::ReductionData)

    Γ0 = canonical_form[:Γ0]
    Γ1 = canonical_form[:Γ1]
    Ψ = canonical_form[:Ψ]
    Π = canonical_form[:Π]
    C = canonical_form[:C]
    n_v = grids[:n_v]
    n_g = grids[:n_g]

    # Perform state space and V reduction, then solve model

    println("performing Krylov reduction on state space...")
    if red_params.reduce_distribution
        Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, basis, basis_inv, n_g_red = krylov_reduction(Γ0, Γ1, Ψ, Π, C, grids, red_params)
    else
        Γ1_kry = Γ1; Γ0_kry = Γ0; Ψ_kry = Ψ; Π_kry = Π; C_kry = C
    end

    println("performing value function reduction via quadratic spline projection...")

    if red_params.reduce_v
        Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline = valuef_reduction(grids, red_params, Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry)
    else
        # Create identity matrix for code reuse below
        from_spline = speye(n_g_red + n_v)
        to_spline = speye(n_g_red + n_v)
        n_splined = n_v
        Γ0_spl = Γ0_kry; Γ1_spl = Γ1_kry; Ψ_spl = Ψ_kry; Π_spl = Π_kry; C_spl = C_kry
    end

    println("solving the model...")
    if red_params.reduce_v && red_params.reduce_distribution
        complex_decomp = false
        inverse_basis = basis_inv*from_spline
    elseif red_params.reduce_distribution
        complex_decomp = false
        inverse_basis = basis_inv
    else
        complex_decomp = true
        inverse_basis = eye(size(T_spl,1))
    end
    T_spl, C_spl, R_spl, ~, ~, ~, ~, eu_spl = gensysct(Γ0_spl, Γ1_spl, C_spl, Ψ_spl, Π_spl,
                                                       complex_decomposition = complex_decomp)

end
