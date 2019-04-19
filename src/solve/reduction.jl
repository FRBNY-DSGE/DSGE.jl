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
function krylov_reduction(m::AbstractCTModel{Float64}, Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                          Ψ::Matrix{Float64}, Π::Matrix{Float64}, C::Array{Float64})

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
    n_state_vars          = n_state_vars - n_state_vars_unreduce

    # Slice Dynamic Equation Γ1 into different parts. See Why Inequality matters paper.
    B_pv = -Γ1[n_total .+ 1:n_vars, n_total .+ 1:n_vars] \ Γ1[n_total .+ 1:n_vars, 1:n_jump_vars]
    B_pg = -Γ1[n_total .+ 1:n_vars, n_total .+ 1:n_vars] \
        Γ1[n_total .+ 1:n_vars, n_jump_vars .+ 1:n_jump_vars + n_state_vars]
    B_pZ = -Γ1[n_total .+ 1:n_vars, n_total .+ 1:n_vars] \
        Γ1[n_total .+ 1:n_vars, n_jump_vars + n_state_vars .+ 1:n_jump_vars + n_state_vars+n_state_vars_unreduce]
    B_gg = Γ1[n_jump_vars .+ 1:n_jump_vars + n_state_vars, n_jump_vars .+ 1:n_jump_vars + n_state_vars]
    B_gv = Γ1[n_jump_vars .+ 1:n_jump_vars + n_state_vars, 1:n_jump_vars]
    B_gp = Γ1[n_jump_vars .+ 1:n_jump_vars + n_state_vars, n_total .+ 1:n_vars]

    # Drop redundant equations in B_pg
    obs        = B_pg
    ~, d0, V_g = svd(obs)
    aux        = d0/d0[1]
    n_Bpg      = Int64(sum(aux .> 10*eps()))
    V_g        = V_g[:, 1:n_Bpg] .* aux[1:n_Bpg]'

    # Compute Krylov subspace
    A(x::Matrix{Float64}) = (F == identity ? 0. : F(B_gv' * x)) + B_gg' * x + B_pg' * (B_gp' * x)
    V_g, ~, ~             = deflated_block_arnoldi(A, V_g, krylov_dim)
    n_state_vars_red      = size(V_g, 2) # number of state variables after reduction

    # Build state space reduction transform
    reduced_basis = spzeros(Float64, n_jump_vars+n_state_vars_red,n_vars)
    reduced_basis[1:n_jump_vars,1:n_jump_vars] = speye(n_jump_vars)
    reduced_basis[n_jump_vars .+ 1:n_jump_vars+n_state_vars_red,n_jump_vars .+ 1:n_jump_vars+n_state_vars] = V_g'
    reduced_basis[n_jump_vars+n_state_vars_red .+ 1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce,n_jump_vars+n_state_vars .+ 1:n_jump_vars+n_state_vars+n_state_vars_unreduce] = eye(n_state_vars_unreduce)

    # Build inverse transform
    inv_reduced_basis = spzeros(n_vars,n_jump_vars+n_state_vars_red)
    inv_reduced_basis[1:n_jump_vars,1:n_jump_vars] = speye(n_jump_vars)
    inv_reduced_basis[n_jump_vars+1:n_jump_vars+n_state_vars,n_jump_vars+1:n_state_vars_red+n_jump_vars] = V_g
    inv_reduced_basis[n_total .+ 1:n_vars,1:n_jump_vars] = B_pv
    inv_reduced_basis[n_total .+ 1:n_vars,n_jump_vars .+ 1:n_jump_vars+n_state_vars_red] = B_pg*V_g
    inv_reduced_basis[n_total .+ 1:n_vars,n_jump_vars+n_state_vars_red .+ 1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = B_pZ
    inv_reduced_basis[n_jump_vars+n_state_vars .+ 1:n_total,n_jump_vars+n_state_vars_red .+ 1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = speye(n_state_vars_unreduce)

    m <= Setting(:n_state_vars_red, n_state_vars_red + n_state_vars_unreduce,
                 "Number of state variables after reduction")

    # change basis
    Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry = change_basis(reduced_basis, inv_reduced_basis, Γ0, Γ1, Ψ, Π, C; ignore_Γ0 = true)

    return Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, reduced_basis, inv_reduced_basis
end

function valuef_reduction(m::AbstractModel,
                          Γ0::Array{Float64,2}, Γ1::Array{Float64,2}, Ψ::Array{Float64,2},
                          Π::Array{Float64,2}, C::Array{Float64})

    n_jump_vars      = Int64(get_setting(m, :n_jump_vars))
    n_state_vars_red = Int64(get_setting(m, :n_state_vars_red))
    spline_grid      = get_setting(m, :spline_grid)
    knots_dict       = Dict{Int64,Vector{Float64}}()

    merge!(knots_dict, get_setting(m, :knots_dict))
    n_prior = Int64(get_setting(m, :n_prior))
    n_post  = Int64(get_setting(m, :n_post))

    # Function calls to create basis reduction
    knots_dim       = Int64(sum(map(i -> length(knots_dict[i]) + 1, 1:length(keys(knots_dict)))))
    spline_grid_dim = Int64(prod(size(spline_grid)))
    from_spline     = spzeros(Float64, spline_grid_dim, knots_dim)
    to_spline       = spzeros(Float64, knots_dim, spline_grid_dim)

    from_spline[:,:], to_spline[:,:] = spline_basis(spline_grid, knots_dict, 2) # create spline basis

    # Extend to other dimensions along which we did not use a spline approximation.
    # In this case, we computed a spline basis for wealth but not for income
    from_spline, to_spline = extend_to_nd(from_spline, to_spline, n_prior, n_post)

    # extra jump variables besides value function
    extra_jv = Int64(n_jump_vars - length(m.endogenous_states[:value_function]))

    if extra_jv > 0
        # Add additional dimensions for any jump variables we did not reduce b/c
        # spline reduction only for value function
        dim1_from, dim2_from = size(from_spline)
        from_spline = sparse([from_spline spzeros(dim1_from, extra_jv)]) # add extra columns
        from_spline = sparse([from_spline; spzeros(extra_jv, dim2_from + extra_jv)])  # add extra zeros
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
    Π_spl = zeros(Float64, Int64(n_states(m)), Int64(n_shocks_expectational(m)))

    Π_spl = to_spline * Π * from_spline[1:n_jump_vars, 1:n_splined]
    m <= Setting(:n_splined, n_splined, "Dimension of jump variables after spline basis reduction")

    return Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline
 end

# Look up the algorithm for this method.
function deflated_block_arnoldi(A::Function, B::Matrix{Float64}, m::Int64)
    F = qr(B)
    Q = Matrix{Float64}(F.Q * Matrix{Float64}(I, size(B)))
    basis = Matrix{Float64}(undef, size(Q,1), 0)
    realsmall = sqrt(eps())

    if m == 1
        basis = Q
    else
        for i in 1:m-1
            # Manual Gram-Schmidt
            basis = hcat(basis, Q)
            aux = A(Q)::Matrix{Float64}
            for j in 1:size(basis,2)
                aux -= basis[:,j] .* basis[:,j]' * aux
            end

            # Check for potential deflation
            Q = Matrix{Float64}(undef, size(Q,1), 0)
            for j in 1:size(aux,2)
                weight = sqrt(sum(aux[:,j].^2))
                if weight > realsmall
                    Q = hcat(Q, aux[:,j] / weight)
                    for k in j+1:size(aux,2)
                        aux[:,k] -= Q[:,end] * (Q[:,end]' * aux[:,k])
                    end
                end
            end

            # More Gram-Schmidt
            for j in 1:size(basis,2)
                Q -= basis[:,j] .* basis[:,j]' * Q
            end
            Q = Q ./ sqrt.(sum(Q.^2, dims=1))
        end
    end
    err = qr(A(Q) - basis * basis' * A(Q)).R::Matrix{Float64}
    return basis, Q, err
end

# Just think through how matrices can be used to do a change of basis
function change_basis(basis, inv_basis, Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                      Ψ::Matrix{Float64}, Π::Matrix{Float64}, C::Array{Float64};
                      ignore_Γ0::Bool = false)
    g1 = basis * Γ1 * inv_basis
    if ignore_Γ0
        g0 = eye(g1)
    else
        g0 = basis * Γ0 * inv_basis
    end
    c    = basis * C
    psi  = basis * Ψ
    Pi   = basis * Π

    return g0, g1, psi, Pi, c
end

# Creates a quadratic polynomial based spline basis reduction
# QUASI - DEPRECATED: translation from SeHyoun's code but not used b/c not general, better
# packages available already in Julia
# Parameters:
#    x = fine grid points
#    knots = coarser grid points onto which we reduce
# Outputs:
#    from_knots = change of basis from spline basis to finer grid points
#    to_knots = change of basis from finer grid points to spline basis
function oneDquad_spline(x::Vector{Float64}, knots::Vector{Float64})
    n_a = length(x)
    n_knots = length(knots)

    # Preallocate auxiliary matrices
    first_interp_mat = zeros(n_a, n_knots + 1)
    aux_mat = zeros(n_a, n_knots)

    # Find to where each points corresponds on knots
    for i = 1:n_a
        loc = sum(knots .<= x[i])
        if loc == n_knots
            loc = n_knots - 1
        end
        first_interp_mat[i, loc] = 1 - (x[i] - knots[loc])^2/(knots[loc + 1] - knots[loc])^2
        first_interp_mat[i, loc + 1] = (x[i] - knots[loc])^2/(knots[loc + 1] - knots[loc])^2
        aux_mat[i, loc] = (x[i] - knots[loc]) - (x[i] - knots[loc])^2/(knots[loc + 1] - knots[loc])
    end

    aux_mat2 = spdiagm(ones(n_knots), 0, n_knots, n_knots) + spdiagm(ones(n_knots-1), 1, n_knots, n_knots)
    aux_mat2[end, end] = 0
    aux_mat2[n_knots, 1] = 1
    aux_mat3 = spdiagm([-2 ./ diff(knots); 0], 0, n_knots, n_knots + 1) +
        spdiagm([2 ./ diff(knots); 1], 1, n_knots, n_knots + 1)

    # Return values
    from_knots = first_interp_mat + aux_mat*(Matrix(aux_mat2)\Matrix(aux_mat3))
    to_knots = from_knots'*from_knots\(from_knots'*eye(n_a))
    return sparse(from_knots), sparse(to_knots)
end

# Creates the spline basis matrix and its inverse, assuming the same degree of splines
# in each direction of approximation
function spline_basis(x::Array{Float64}, knots::Dict{Int64,Vector{Float64}}, degree::Int64)

    # Create spline approximations along each dimension
    approx_dim = length(size(x))
    basis = Basis(map(i -> Basis(SplineParams(knots[i], 0, degree)), 1:approx_dim)...)

    # Create spline basis matrix and its inverse. The spline basis matrix transforms
    # the coefficient vector of knot points to the coefficient vector of the provided
    # state space array x, i.e. S*knots = x, where S is the spline basis matrix.
    from_knots = BasisMatrix(basis, Expanded(), x, 0).vals[1]
    to_knots = Matrix{Float64}(from_knots'*from_knots) \ Matrix{Float64}(from_knots')
    return from_knots, sparse(to_knots)
end

# Creates the spline basis matrix and its inverse, assuming a vector specifying the desired degrees.
function spline_basis(x::Array{Float64}, knots::Dict{Int64,Vector{Float64}}, degree::Vector{Int64})
    # Create spline approximations along each dimension
    approx_dim = length(size(x))
    basis = Basis(SplineParams(knots[1], 0, degree))
    for i = 2:approx_dim
        basis = Basis(basis, SplineParams(knots[i], 0, degree)) # UNSATISFIED WITH THIS; FASTER WAY?
    end

    # Create spline basis matrix and its inverse. The spline basis matrix transforms
    # the coefficient vector of knot points to the coefficient vector of the provided
    # state space array x, i.e. S*knots = x, where S is the spline basis matrix.
    from_knots = BasisMatrix(basis, Expanded(), x, 0).vals[1]
    to_knots = Matrix(from_knots'*from_knots) \ Matrix(from_knots')
    return from_knots, sparse(to_knots)
end

# Parameters:
#    from_small = projection to approximation
#    to_small = approximation to projection
#    n_prior = number of grid points stacked prior to reduction dimension
#    n_post = number of grid points stacked after reduction dimension
#
# Outputs:
#    from_approx = basis change from small points to large
#    to_approx = basis change from large number of points to small
function extend_to_nd(from_small::SparseMatrixCSC{Float64,Int64}, to_small::SparseMatrixCSC{Float64,Int64},
                      n_prior::Int64, n_post::Int64)
    from_approx = kron(kron(SparseMatrixCSC{Float64}(I, n_post, n_post), from_small),
                       SparseMatrixCSC{Float64}(I, n_prior, n_prior))
    to_approx = kron(kron(SparseMatrixCSC{Float64}(I, n_post, n_post), to_small),
                     SparseMatrixCSC{Float64}(I, n_prior, n_prior))
    return from_approx, to_approx
end

# Create transformations when we only transform subset of the points.
#
# Arguments
#
# from_small: (nxk) matrix transformation from approximated k points into full grids
# to_small: (kxn) matrix transformation into approximation with k points
# n_pre: number of points preceding the transformation
# n_post: number of points following the transformation
#
# Outputs:
# from_approx: transformation from smaller approximation to larger grid points
# to_approx: reverse of from_approx
function projection_for_subset(from_small::SparseMatrixCSC{Float64,Int64},
                               to_small::SparseMatrixCSC{Float64,Int64},
                               n_pre::Int64, n_post::Int64)

    n_full, n_red = size(from_small)
    I, J, V = SparseArrays.spdiagm_internal(0 => ones(n_red+n_pre)); from_approx = sparse(I, J, V, n_full + n_pre, n_pre + n_red)
    I, J, V = SparseArrays.spdiagm_internal(0 => ones(n_red+n_pre)); to_approx = sparse(I, J, V, n_pre + n_red, n_full + n_pre)
    from_approx[1+n_pre:n_pre+n_full, n_pre+1:n_pre+n_red] = from_small
    to_approx[n_pre+1:n_pre+n_red, n_pre+1:n_pre+n_full] = to_small

    # Expand matrices and add needed values
    dim1_from_approx, dim2_from_approx = size(from_approx)
    from_approx = [from_approx spzeros(dim1_from_approx, n_post)]
    from_approx = [from_approx; spzeros(n_post, dim2_from_approx + n_post)]
    to_approx = [to_approx spzeros(dim2_from_approx, n_post)]
    to_approx = [to_approx; spzeros(n_post, dim1_from_approx + n_post)]
    from_approx[dim1_from_approx+1:end, dim2_from_approx+1:end] = speye(n_post)
    to_approx[dim2_from_approx+1:end, dim1_from_approx+1:end] = speye(n_post)

    return from_approx, to_approx
end

# Cleans Γ0 matrix to be identity by solving out static conditions
function solve_static_conditions(Γ0::Matrix{Float64},
                                 Γ1::Matrix{Float64},
                                 Ψ::Matrix{Float64},
                                 Π::Matrix{Float64},
                                 C::Array{Float64})

    redundant = maximum(abs.([Γ0 Ψ]), 2) .== 0  # Find rows of both Γ0 & Ψ with only zeros
    redundant_list = Array{Int64,1}(0)
    for row in 1:length(redundant)
        if redundant[row]
            push!(redundant_list, row)
        end
    end
    inv_state_red = nullspace(Γ1[redundant_list,:]) # Compute orthonormal basis of null space of Γ1
    state_red = inv_state_red'

    g0  = state_red * Γ0 * inv_state_red # adjust Γ0 to new basis
    g1  = state_red * Γ1 * inv_state_red # Zero out redundant rows/columns in Γ1
    g1  = g0 \ g1                        # Force Γ1 and other matrices
    Psi = g0 \ (state_red * Ψ)           # to be the ones corresponding to when
    Pi  = g0 \ (state_red * Π)           # Γ0 is the identity
    c   = g0 \ (state_red * C)

    return g0, g1, c, Pi, Psi, state_red, inv_state_red
end

# When we have sparse matrices
function solve_static_conditions(Γ0::SparseMatrixCSC{Float64, Int64},
                                 Γ1::SparseMatrixCSC{Float64, Int64},
                                 Ψ::SparseMatrixCSC{Float64, Int64},
                                 Π::SparseMatrixCSC{Float64, Int64},
                                 C::SparseMatrixCSC{Float64, Int64})

    redundant = maximum(abs.([Γ0 Ψ]), 2) .== 0 # R: I removed the comma here, given the spirit of what
                                               # it appears we're trying to do.
    redundant_list = Array{Int64,1}(0)
    for row in 1:length(redundant)
        if redundant[row]
            push!(redundant_list, row)
        end
    end
    inv_state_red = nullspace(full(Γ1[redundant_list,:]))
    state_red = inv_state_red'

    g0 = state_red * Matrix(Γ0) * inv_state_red
    g1 = state_red * Matrix(Γ1) * inv_state_red
    g1 = g0 \ g1
    Psi = g0 \ (state_red * Matrix(Ψ))
    Pi = g0 \ (state_red * Matrix(Π))
    c = g0 \ (state_red * Vector(C))

    return g0, g1, c, Pi, Psi, state_red, inv_state_red
end

# When we have sparse matrices and a vector for C
function solve_static_conditions(Γ0::SparseMatrixCSC{Float64, Int64},
                                 Γ1::SparseMatrixCSC{Float64, Int64},
                                 Ψ::SparseMatrixCSC{Float64, Int64},
                                 Π::SparseMatrixCSC{Float64, Int64},
                                 C::SparseVector{Float64, Int64})

    redundant = maximum(abs.([Γ0, Ψ]), 2) .== 0
    inv_state_red = null(full(Γ1(redundant,:)))
    state_red = inv_state_red'

    g0 = state_red * Matrix(Γ0) * inv_state_red
    g1 = state_red * Matrix(Γ1) * inv_state_red
    g1 = g0 \ g1
    Psi = g0 \ (state_red * Matrix(Ψ))
    Pi = g0 \ (state_red * Matrix(Π))
    c = g0 \ (state_red * Vector(C))

    return g0, g1, c, Pi, Psi, state_red, inv_state_red
end
