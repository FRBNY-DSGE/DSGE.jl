"""
```
krylov_reduction(Γ0, Γ1, Ψ, Π, C; F = x -> x, n_Z = 0)
valuef_reduction(grids, params, n_g_red, Γ0, Γ1, Ψ, Π, C, n_v)
```

Perform Krylov reduction on the state space system defined by Γ0 and Γ1
and/or value function reduction using a spline approximation across some
subset of the state space.

Returns the new matrices, the change of basis matrix and its inverse,
and the new size of the state space.

"""
function krylov_reduction(m::AbstractCTModel, Γ0::SparseMatrixCSC{Float64,Int64},
                          Γ1::SparseMatrixCSC{Float64,Int64},
                          Ψ::SparseMatrixCSC{Float64,Int64}, Π::SparseMatrixCSC{Float64,Int64},
                          C::SparseMatrixCSC{Float64,Int64})

    # Get Krylov subspace function
    if !in(:F, keys(m.settings))
        m <= Setting(:F, identity, "Function applied during Krylov reduction")
    end
    F = get_setting(m, :F)

    # Grab Dimensions
    n_state_vars_unreduce = Int64(get_setting(m, :n_state_vars_unreduce))
    n_state_vars          = Int64(get_setting(m, :n_state_vars))
    n_jump_vars           = Int64(get_setting(m, :n_jump_vars))
    krylov_dim            = Int64(get_setting(m, :krylov_dim))
    n_total               = n_jump_vars + n_state_vars
    n_vars                = n_states(m)
    n_state_vars         -= n_state_vars_unreduce

    # Slice Γ1 into different parts. See Why Inequality matters paper.
    Γ1_full = Array{Float64}(Γ1)
    B_pv    = -Γ1_full[n_total .+ 1:n_vars, n_total + 1:n_vars] \ Γ1_full[n_total .+ 1:n_vars, 1:n_jump_vars]
    B_pg    = -Γ1_full[n_total .+ 1:n_vars, n_total .+ 1:n_vars] \
           Γ1_full[n_total .+ 1:n_vars, n_jump_vars .+ 1:n_jump_vars + n_state_vars]
    B_pZ    = -Γ1_full[n_total .+ 1:n_vars, n_total .+ 1:n_vars] \
        Γ1_full[n_total+1:n_vars, n_jump_vars+n_state_vars .+ 1:n_jump_vars+n_state_vars+n_state_vars_unreduce]
    B_gg    = Γ1_full[n_jump_vars .+ 1:n_jump_vars + n_state_vars, n_jump_vars .+ 1:n_jump_vars + n_state_vars]
    B_gv    = Γ1_full[n_jump_vars .+ 1:n_jump_vars + n_state_vars, 1:n_jump_vars]
    B_gp    = Γ1_full[n_jump_vars .+ 1:n_jump_vars .+ n_state_vars, n_total .+ 1:n_vars]
    # Drop redundant equations
    obs        = B_pg
    ~, d0, V_g = svd(obs)
    aux        = d0/d0[1]
    n_Bpg      = Int64(sum(aux .> 10*eps()))
    V_g        = V_g[:, 1:n_Bpg] .* aux[1:n_Bpg]'
    # Compute Krylov subspace
    A(x::Matrix{Float64}) = (F == identity ? 0. : F(B_gv' * x)) .+ B_gg' * x + B_pg' * (B_gp' * x)
    V_g, ~, ~             = deflated_block_arnoldi(A, V_g, krylov_dim)
    n_state_vars_red      = size(V_g, 2) # number of state variables after reduction
    # Build state space reduction transform
    reduced_basis = spzeros(Float64, n_jump_vars+n_state_vars_red,n_vars)
    reduced_basis[1:n_jump_vars,1:n_jump_vars] = SparseMatrixCSC{Float64}(I, n_jump_vars, n_jump_vars)
    reduced_basis[n_jump_vars .+ 1:n_jump_vars+n_state_vars_red,n_jump_vars .+ 1:n_jump_vars+n_state_vars] = V_g'
    reduced_basis[n_jump_vars + n_state_vars_red .+ 1:n_jump_vars + n_state_vars_red + n_state_vars_unreduce,
                  n_jump_vars + n_state_vars .+ 1:n_jump_vars + n_state_vars + n_state_vars_unreduce] =
                      Matrix{Float64}(I, n_state_vars_unreduce,
                                      n_state_vars_unreduce)

    # Build inverse transform
    inv_reduced_basis                              = spzeros(Float64, n_vars, n_jump_vars + n_state_vars_red)
    inv_reduced_basis[1:n_jump_vars,1:n_jump_vars] = SparseMatrixCSC(I, n_jump_vars, n_jump_vars)
    inv_reduced_basis[n_jump_vars+1:n_jump_vars+n_state_vars,n_jump_vars+1:n_state_vars_red+n_jump_vars] = V_g
    inv_reduced_basis[n_total+1:n_vars,1:n_jump_vars] = B_pv
    inv_reduced_basis[n_total+1:n_vars,n_jump_vars+1:n_jump_vars .+ n_state_vars_red] = B_pg * V_g
    inv_reduced_basis[n_total+1:n_vars,n_jump_vars+n_state_vars_red .+ 1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = B_pZ
    inv_reduced_basis[n_jump_vars+n_state_vars .+ 1:n_total,n_jump_vars+n_state_vars_red .+ 1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = SparseMatrixCSC{Float64}(I, n_state_vars_unreduce, n_state_vars_unreduce)
    m <= Setting(:n_state_vars_red, n_state_vars_red + n_state_vars_unreduce,
                 "Number of state variables after reduction")

    # change basis
    Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry =
        change_basis(reduced_basis, inv_reduced_basis, Γ0, Γ1, Ψ, Π, C; ignore_Γ0 = true)

    return Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, reduced_basis, inv_reduced_basis
end

function valuef_reduction(m::AbstractCTModel,
                          Γ0::SparseMatrixCSC{Float64,Int64}, Γ1::SparseMatrixCSC{Float64,Int64},
                          Ψ::SparseMatrixCSC{Float64,Int64},
                          Π::SparseMatrixCSC{Float64,Int64}, C::SparseMatrixCSC{Float64,Int64})

    n_jump_vars      = Int64(get_setting(m, :n_jump_vars))
    n_state_vars_red = Int64(get_setting(m, :n_state_vars_red))
    spline_grid      = get_setting(m, :spline_grid)
    knots_dict       = get_setting(m, :knots_dict)
    n_prior          = Int64(get_setting(m, :n_prior))
    n_post           = Int64(get_setting(m, :n_post))

    # Function calls to create basis reduction
    knots_dim       = sum(map(i -> length(knots_dict[i]) + 1, 1:length(keys(knots_dict))))
    spline_grid_dim = prod(size(spline_grid))
    from_spline     = spzeros(Float64, spline_grid_dim, knots_dim)
    to_spline       = spzeros(Float64, knots_dim, spline_grid_dim)

    # create spline basis
    from_spline[:,:], to_spline[:,:] = spline_basis(spline_grid, knots_dict, 2)
    # extend to other dimensions along which we did not use a spline approximation.
    # In this case, we computed a spline basis for wealth but not for income
    from_spline, to_spline = extend_to_nd(from_spline, to_spline, n_prior, n_post)
    # extra jump variables besides value function
    extra_jv = Int64(n_jump_vars - length(m.endogenous_states[:value_function]))

    if extra_jv > 0
        # Add additional dimensions for any jump variables we did not reduce b/c
        # spline reduction only for value function
        dim1_from, dim2_from = size(from_spline)
        from_spline = sparse([from_spline spzeros(dim1_from, extra_jv)]) # add extra columns
        from_spline = sparse([from_spline; spzeros(extra_jv, dim2_from + extra_jv)]) # add extra zeros
        to_spline = sparse([to_spline spzeros(dim2_from, extra_jv)])     # same as above
        to_spline = sparse([to_spline; spzeros(1, dim1_from + extra_jv)])
        for i = 0:(extra_jv - 1)
            from_spline[end - i, end - i] = 1 # add 1 along diagonal of added dimensions
            to_spline[end - i, end - i] = 1
        end
    end
    n_splined = Int64(size(from_spline, 2))

    # Create projection matrix that projects value function onto spline basis
    from_spline, to_spline = projection_for_subset(from_spline, to_spline, 0, n_state_vars_red)

    # Reduce the decision vector
    Γ0_spl, Γ1_spl, Ψ_spl, ~, C_spl = change_basis(to_spline, from_spline, Γ0, Γ1, Ψ, Π, C; ignore_Γ0 = true)

    Π_spl = to_spline * Π * from_spline[1:n_jump_vars, 1:n_splined]
    m  <= Setting(:n_splined, n_splined, "Dimension of jump variables after spline basis reduction")

    return Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline
 end

# Just think through how matrices can be used to do a change of basis
function change_basis(basis, inv_basis, Γ0::SparseMatrixCSC{Float64,Int64},
                      Γ1::SparseMatrixCSC{Float64,Int64},
                      Ψ::SparseMatrixCSC{Float64,Int64}, Π::SparseMatrixCSC{Float64,Int64},
                      C::SparseMatrixCSC{Float64,Int64};
                      ignore_Γ0::Bool = false)

    g1 = basis * Γ1 * inv_basis

    if ignore_Γ0
        g0 = SparseMatrixCSC{Float64}(I, size(g1))
    else
        g0 = basis * Γ0 * inv_basis
    end
    c    = basis * C
    psi  = basis * Ψ
    Pi   = basis * Π

    return g0, g1, psi, Pi, c
end

# Cleans Γ0 matrix to be identity by solving out static conditions
function solve_static_conditions(Γ0::SparseMatrixCSC{Float64,Int64},
                                 Γ1::SparseMatrixCSC{Float64,Int64},
                                 Ψ::SparseMatrixCSC{Float64,Int64},
                                 Π::SparseMatrixCSC{Float64,Int64}, C::Array{Float64})

    redundant     = maximum(abs.([Γ0, Ψ]), 2) .== 0  # Find rows of both Γ0 & Ψ with only zeros
    inv_state_red = null(Γ1(redundant,:)) # Compute orthonormal basis of null space of Γ1
    state_red     = inv_state_red'

    g0  = state_red * Γ0 * inv_state_red # adjust Γ0 to new basis
    g1  = state_red * Γ1 * inv_state_red # Zero out redundant rows/columns in Γ1
    g1  = g0 \ g1                        # Force Γ1 and other matrices
    Psi = g0 \ (state_red * Ψ)           # to be the ones corresponding to when
    Pi  = g0 \ (state_red * Π)           # Γ0 is the identity
    c   = g0 \ (state_red * C)

    return g0, g1, c, Pi, Psi, state_red, inv_state_red
end

# When we have sparse matrices
#=
function solve_static_conditions(Γ0::SparseMatrixCSC{Float64, Int64}, Γ1::SparseMatrixCSC{Float64, Int64},
                  Ψ::SparseMatrixCSC{Float64, Int64},
                  Π::SparseMatrixCSC{Float64, Int64}, C::SparseMatrixCSC{Float64, Int64})

    redundant     = maximum(abs.([Γ0, Ψ]), 2) .== 0
    inv_state_red = null(full(Γ1(redundant,:)))
    state_red     = inv_state_red'

    g0  = state_red * full(Γ0) * inv_state_red
    g1  = state_red * full(Γ1) * inv_state_red
    g1  = g0 \ g1
    Psi = g0 \ (state_red * full(Ψ))
    Pi  = g0 \ (state_red * full(Π))
    c   = g0 \ (state_red * full(C))

    return g0, g1, c, Pi, Psi, state_red, inv_state_red
end
=#
