"""
```
solvect(m::AbstractCTModel; reduction_settings_check::Bool = false)
```

Driver to compute the model solution and augment transition matrices.

### Inputs

- `m`: the model object

### Outputs
 - TTT, RRR, and CCC matrices of the state transition equation:
    ```
    S_t = TTT*S_{t-1} + RRR*ϵ_t + CCC
    ```
"""
function solve(m::AbstractCTModel; sparse_mat::Bool = true, check_Γ0::Bool = true,
               complex_decomp::Bool = false, alt_policy::Bool = false)

    # Get equilibrium condition matrices
    Γ0, Γ1, Ψ, Π, C = eqcond(m)

    # Check that Γ0 is approximately the identity, otherwise clean it by solving static conditions
    if check_Γ0
        try
            n_total = get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars)
            @assert Matrix{Float64}(I, n_total, n_total) ≈ Γ0[1:n_total, 1:n_total]
        catch
            Γ0, Γ1, Ψ, Π, C, basis_redundant, inv_basis_redundant = solve_static_conditions(Γ0, Γ1, Ψ, Π, C)
        end
    end
    if !(@isdefined inv_basis_redundant)
        inv_basis_redundant = SparseMatrixCSC{Float64}(I, size(Γ0,1), size(Γ0,1))
        basis_redundant = SparseMatrixCSC{Float64}(I, size(Γ0,1), size(Γ0,1))
    end

    if sparse_mat
        Γ0 = sparse(Γ0); Γ1 = sparse(Γ1); Ψ = sparse(Ψ); Π = sparse(Π)
        C = sparse(reshape(C, length(C), 1))
    end

    # krylov reduction
    if get_setting(m, :reduce_state_vars)
        Γ0, Γ1, Ψ, Π, C, basis_kry, inv_basis_kry = krylov_reduction(m, Γ0, Γ1, Ψ, Π, C)
    end

    # value function reduction via spline projection
    if get_setting(m, :reduce_v)
        Γ0, Γ1, Ψ, Π, C, basis_spl, inv_basis_spl = valuef_reduction(m, Γ0, Γ1, Ψ, Π, C)
    end

    # Compute inverse basis for Z matrix and IRFs transformation
    if get_setting(m, :reduce_v) && get_setting(m, :reduce_state_vars)
        inverse_basis = inv_basis_redundant * inv_basis_kry * inv_basis_spl
        basis = basis_spl * basis_kry * basis_redundant
    elseif get_setting(m, :reduce_state_vars)
        inverse_basis = inv_basis_redundant * inv_basis_kry
        basis = basis_kry * basis_redundant
    else
        inverse_basis = inv_basis_redundant
        basis = basis_redundant
    end

    # Solve LRE model
    TTT_gensys, CCC_gensys, RRR_gensys, ~, ~, eu = gensysct!(Matrix{Float64}(Γ1), Matrix{Float64}(C),
                                                             Matrix{Float64}(Ψ), Matrix{Float64}(Π),
                                                             complex_decomposition = complex_decomp)

    # Check for LAPACK exception, existence and uniqueness
    if eu[1] != 1 || eu[2] != 1
        throw(GensysError())
    end

    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = vec(real(CCC_gensys))

    # If you wanted to augment states, you would do it here

    return TTT_gensys, RRR_gensys, CCC_gensys, inverse_basis, basis
end
