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
function krylov_reduction(Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                          Ψ::Matrix{Float64}, Π::Matrix{Float64}, C::Array{Float64},
                          grid::Dict{Symbol, Any}, red_params::ReductionData;
                          F::Function = identity, n_state_vars_unreduce = 0)

    # Grab Dimensions
    n_state_vars = grids[:n_state_vars]
    n_jump_vars = grids[:n_jump_vars]
    n_p = grids[:n_static_conditions]
    krylov_dim = red_params.krylov_dim
    n_total = n_jump_vars + n_state_vars
    n_state_vars = n_state_vars - n_state_vars_unreduce

    if !(eye(n_total) ≈ Γ0[1:n_total, 1:n_total])
        error("Γ0 is not normalized")
    end

    # Slice Γ1 into different parts
    B_pv = -Γ1[n_total+1:n_total+n_p,n_total+1:n_total+n_p]\Γ1[n_total+1:n_total+n_p,1:n_jump_vars]
    B_pg = -Γ1[n_total+1:n_total+n_p,n_total+1:n_total+n_p]\Γ1[n_total+1:n_total+n_p,n_jump_vars+1:n_jump_vars+n_state_vars]
    B_pZ = -Γ1[n_total+1:n_total+n_p,n_total+1:n_total+n_p]\Γ1[n_total+1:n_total+n_p,n_jump_vars+n_state_vars+1:n_jump_vars+n_state_vars+n_state_vars_unreduce]
    B_gg = Γ1[n_jump_vars+1:n_jump_vars+n_state_vars,n_jump_vars+1:n_jump_vars+n_state_vars]
    B_gv = Γ1[n_jump_vars+1:n_jump_vars+n_state_vars,1:n_jump_vars]
    B_gp = Γ1[n_jump_vars+1:n_jump_vars+n_state_vars,n_total+1:n_total+n_p]

    # Drop redundant equations
    obs = B_pg
    ~, d0, V_g = svd(obs)
    aux = d0/d0[1]
    n_Bpg = sum(aux .> 10*eps())
    V_g = V_g[:, 1:n_Bpg] .* aux[1:n_Bpg]'

    # Compute Krylov subspace
    A(x::Matrix{Float64}) = (F == identity ? 0. : F(B_gv' * x)) + B_gg' * x + B_pg' * (B_gp' * x)
    V_g, ~, ~ = deflated_block_arnoldi(A, V_g, krylov_dim)
    n_state_vars_red = size(V_g, 2)

    # Build state space reduction transform
    reduced_basis = spzeros(n_jump_vars+n_state_vars_red,n_total+n_p)
    reduced_basis[1:n_jump_vars,1:n_jump_vars] = speye(n_jump_vars)
    reduced_basis[n_jump_vars+1:n_jump_vars+n_state_vars_red,n_jump_vars+1:n_jump_vars+n_state_vars] = V_g'
    reduced_basis[n_jump_vars+n_state_vars_red+1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce,n_jump_vars+n_state_vars+1:n_jump_vars+n_state_vars+n_state_vars_unreduce] = eye(n_state_vars_unreduce)

    # Build inverse transform
    inv_reduced_basis = spzeros(n_total+n_p,n_jump_vars+n_state_vars_red)
    inv_reduced_basis[1:n_jump_vars,1:n_jump_vars] = speye(n_jump_vars)
    inv_reduced_basis[n_jump_vars+1:n_jump_vars+n_state_vars,n_jump_vars+1:n_state_vars_red+n_jump_vars] = V_g
    inv_reduced_basis[n_total+1:n_total+n_p,1:n_jump_vars] = B_pv
    inv_reduced_basis[n_total+1:n_total+n_p,n_jump_vars+1:n_jump_vars+n_state_vars_red] = B_pg*V_g
    inv_reduced_basis[n_total+1:n_total+n_p,n_jump_vars+n_state_vars_red+1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = B_pZ
    inv_reduced_basis[n_jump_vars+n_state_vars+1:n_total,n_jump_vars+n_state_vars_red+1:n_jump_vars+n_state_vars_red+n_state_vars_unreduce] = speye(n_state_vars_unreduce)

    grids[:n_state_vars_red] = n_state_vars_red + n_state_vars_unreduce

    # change basis
    Γ0, Γ1, Ψ, Π, C = change_basis(reduced_basis, inv_reduced_basis, Γ0, Γ1, Ψ, Π, C)

    return Γ0, Γ1, Ψ, Π, C, reduced_basis, inv_reduced_basis
end

function valuef_reduction(grids::Dict{Symbol,Any}, red_params::ReductionData,
                          Γ0::Array{Float64,2}, Γ1::Array{Float64,2}, Ψ::Array{Float64,2},
                          Π::Array{Float64,2}, C::Array{Float64,2})
    n_jump_vars = grids[:n_jump_vars]
    n_state_vars_red = grids[:n_state_vars_red]
    ss_array = red_params.ss_array
    knots_dict = red_params.knots_dict
    n_prior = red_params.n_prior
    n_post = red_params.n_post

    # Function calls to create basis reduction
    from_spline, to_spline = spline_basis(ss_array, knots_dict, 2)
    from_spline, to_spline = extend_to_nd(from_spline, to_spline, n_prior, n_post)
    extra_jv = grids[:n_jump_vars] - length(grids[:a])*length(grids[:z]) # extra jump variables besides state space grid
    if extra_jv > 0
        dim1_from, dim2_from = size(from_spline)
        from_spline = [from_spline spzeros(dim1_from, extra_jv)]
        from_spline = [from_spline; spzeros(extra_jv, dim2_from + extra_jv)]
        to_spline = [to_spline spzeros(dim2_from, extra_jv)]
        to_spline = [to_spline; spzeros(1, dim1_from + extra_jv)]
        for i = 0:(extra_jv - 1)
            from_spline[end - i, end - i] = 1
            to_spline[end - i, end - i] = 1
        end
    end
    n_splined = size(from_spline, 2)
    from_spline, to_spline = projection_for_subset(from_spline, to_spline, 0, n_state_vars_red)

    # Reduce the decision vector
    Γ0_spl, Γ1_spl, Ψ_spl, ~, C_spl = change_basis(to_spline, from_spline, Γ0, Γ1, Ψ, Π, C)
    Π_spl = to_spline * Π * from_spline[1:n_jump_vars, 1:n_splined]
    grids[:n_splined] = n_splined

    return Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline
 end

function deflated_block_arnoldi(A::Function, B::Matrix{Float64}, m::Int64)

    Q, ~ = qr(B)
    basis = Matrix{Float64}(size(Q,1), 0)
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
            Q = Matrix{Float64}(size(Q,1), 0)
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
            Q = Q ./ sqrt.(sum(Q.^2,1))
        end
    end

    err = qr(A(Q) - basis * basis' * A(Q))[2]::Matrix{Float64}

    return basis, Q, err
end

function change_basis(basis, inv_basis, Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                      Ψ::Matrix{Float64}, Π::Matrix{Float64}, C::Array{Float64})
    Γ1 = basis * Γ1 * inv_basis
    Γ0 = basis * Γ0 * inv_basis
    C  = basis * C
    Ψ  = basis * Ψ
    Π  = basis * Π

    return Γ0, Γ1, Ψ, Π, C
end

# Creates a quadratic polynomial based spline basis reduction
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
    aux_mat3 = spdiagm([-2./diff(knots); 0], 0, n_knots, n_knots + 1) +
        spdiagm([2./diff(knots); 1], 1, n_knots, n_knots + 1)

    # Return values
    from_knots = first_interp_mat + aux_mat*(full(aux_mat2)\full(aux_mat3))
    to_knots = from_knots'*from_knots\(from_knots'*eye(n_a))
    return sparse(from_knots), sparse(to_knots)
end

# Parameters:
#    from_small = projection to approximation
#    to_small = approximation to projection
#    n_prior = number of grid points stacked prior to reduction dimension
#    n_post = number of gird points stacked after reduction dimension
#
# Outputs:
#    from_approx = basis change from small points to large
#    to_approx = basis change from large number of points to small
function extend_to_nd(from_small::SparseMatrixCSC{Float64,Int64}, to_small::SparseMatrixCSC{Float64,Int64},
                      n_prior::Int64, n_post::Int64)
    from_approx = kron(kron(speye(n_post), from_small), speye(n_prior))
    to_approx = kron(kron(speye(n_post), to_small), speye(n_prior))
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
    from_approx = spdiagm(ones(n_red+n_pre), 0, n_full + n_pre, n_pre + n_red)
    to_approx = spdiagm(ones(n_red+n_pre), 0, n_pre + n_red, n_full + n_pre)
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
    to_knots = full(from_knots'*from_knots) \ full(from_knots')
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
    to_knots = full(from_knots'*from_knots) \ full(from_knots')
    return from_knots, sparse(to_knots)
end

