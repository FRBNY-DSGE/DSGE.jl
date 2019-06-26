# This file contains various methods for solving the
# Kolmogorov-Forward equation in the steady state
# computation for HANK models:
# 0 = A' * g
# They can be generally classified as
# methods for well-conditioned A matrices and ill-conditioned A matrices.
# All of these methods return the stacked distribution.

function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, partitions::Vector{Float64},
                   dimensions::Vector{Int64})

    if size(partitions, 1) != 1 || size(dimensions, 1) != 1
        error("Need to enter partitions and/or dimensions as a 1 x n array or as a vector")
    end

    AT       = A'
    dim_prod = prod(dimensions) # Need total size of dimensions

    # Fixing one value to ensure matrix is not singular, may require further generalization later
    i_fix    = 1
    b        = zeros(dim_prod)
    b[i_fix] = .1
    for j = 1:dim_prod
        AT[i_fix,j] = 0
    end
    AT[i_fix, i_fix] = 1

    # Solve linear system for distribution
    gg    = AT \ b
    g_sum = gg'*ones(dim_prod, 1) * prod(partitions)
    gg    = gg ./ g_sum
    return vec(gg)
end

# Iterative method for ill-conditioned A
# g0: initial distribution guess provided by user
# integrate_mat: weights distribution points by some multiple of their Lebesgue measure across wealth
function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Matrix{Float64},
                   weight_mat::SparseMatrixCSC{Float64, Int64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) != size(A)
        error("Dimension of weight matrix is incorrect.")
    end
    dim_size = size(A, 1)
    gg = vec(g0) # stack distribution matrix into vector

    # Solve linear system
    for ikfe = 1:maxit_kfe
        gg_tilde = weight_mat * gg # weight distribution points by their measure across wealth

        # 0 = A'g -> g = A'g + g -> (I - A')g = g, but possible singularity so perturb from zero
        # by adding g/Δ_KFE, where Δ_KFE is large => g/Δ_KFE ≈ 0. Thus,
        # g/Δ_KFE = A'*g + g/Δ_KFE => g = Δ_KFE A'*g + g
        # => (I - Δ_KFE A')g = g; could also use (I + Δ_KFE A')g = g
        gg1_tilde = (speye(ComplexF64, dim_size) - Δ_kfe * A') \ gg_tilde
        gg1_tilde = gg1_tilde ./ sum(gg1_tilde) # normalize
        gg1       = weight_mat \ gg1_tilde # undo integration over wealth

        # Check iteration for convergence
        err_kfe = maximum(abs.(gg1 - gg))
        if err_kfe < tol_kfe
            break
        end
        gg = gg1
    end

    return vec(gg)
end

# Previous function but for complex g0
function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Matrix{ComplexF64},
                   weight_mat::SparseMatrixCSC{Float64, Int64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) != size(A)
        error("Dimension of weight matrix is incorrect.")
    end
    dim_size = size(A, 1)
    gg = vec(g0) # stack distribution matrix into vector

    # Solve linear system
    for ikfe = 1:maxit_kfe
        gg_tilde = weight_mat * gg # weight distribution points by their measure across wealth

        # 0 = A'g -> g = A'g + g -> (I - A')g = g, but to avoid issues with small numbers
        # we use a large scalar Delta_KFE such that
        # 0 = A' Delta_KFE g -> (I - Delta_KFE A') g = g
        gg1_tilde = (speye(ComplexF64, dim_size) - Δ_kfe * A') \ gg_tilde # Left divide
        gg1_tilde = gg1_tilde ./ sum(gg1_tilde) # normalize to one to get pdf
        gg1 = weight_mat \ gg1_tilde # undo integration over wealth

        # Check iteration for convergence
        err_kfe = maximum(abs.(gg1 - gg))
        if err_kfe < tol_kfe
            break
        end

        gg = gg1
    end

    return vec(gg)
end

# Iterative methods when weight_mat is (1) not sparse or
# (2) a dense matrix that integrates out non-wealth variables (rather than a diagonal one),
function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Vector{Float64},
                   weight_mat::Matrix{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) == size(A)
        return solve_kfe(A, g0, sparse(weight_mat); maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    else
        weight_mat = spdiagm(vec(weight_mat), 0)
        return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    end
end

# Iterative method when weight_mat is a vector of differences of the wealth grid
# e.g. if x = [0 .05, .1, .15], then weight_mat = [.05, .05, .05]
function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Vector{Float64},
                   grid::Vector{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    reps = size(A, 1) / length(grid)
    if Float64(reps) != Float64(Int64(reps))
        error("Dimension of A does not match dimensions of the grid.")
    end
    diagonal = zeros(size(A, 1))
    for i = 1:reps
        reps[1 + length(grid) * (i - 1): length(grid) * i] = grid
    end
    weight_mat = spdiagm(diagonal, 0)
    return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
end

# Allows for complex g0
function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Vector{ComplexF64},
                   weight_mat::Matrix{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) == size(A)
        return solve_kfe(A, g0, sparse(weight_mat); maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    else
        weight_mat = spdiagm(vec(weight_mat), 0)
        return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    end
end

function solve_kfe(A::SparseMatrixCSC{ComplexF64, Int64}, g0::Vector{ComplexF64},
                   grid::Vector{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    reps = size(A, 1) / length(grid)
    if Float64(reps) != Float64(Int64(reps))
        error("Dimension of A does not match dimensions of the grid.")
    end
    diagonal = zeros(size(A, 1))
    for i = 1:reps
        reps[1 + length(grid) * (i - 1): length(grid) * i] = grid
    end
    weight_mat = spdiagm(diagonal, 0)
    return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
end

# Repeats previous code when entries of A are floats.
# Exact method for well-conditioned A
function solve_kfe(A::SparseMatrixCSC{Float64, Int64}, partitions::Vector{Float64},
                   dimensions::Vector{Int64})
    AT = A'
    dim_prod = prod(dimensions)
    # Fixing one value to ensure matrix is not singular
    i_fix = 1
    b = zeros(dim_prod, 1)
    b[i_fix] = .1
    for j = 1:dim_prod
        AT[i_fix,j] = 0
    end
    AT[i_fix, i_fix] = 1

    # Solve linear system for distribution
    gg = AT \ b
    g_sum = gg'*ones(dim_prod, 1) * prod(partitions)
    gg = gg ./ g_sum
    return vec(gg)
end
function solve_kfe(A::SparseMatrixCSC{Float64, Int64}, partitions::Matrix{Float64},
                   dimensions::Matrix{Int64})
    if size(partitions, 1) != 1 || size(dimensions, 1) != 1
        error("Need to enter partitions and/or dimensions as a 1 x n array or as a vector")
    end
    AT = A'
    dim_prod = prod(dimensions)
    # Fixing one value to ensure matrix is not singular
    i_fix = 1
    b = zeros(dim_prod, 1)
    b[i_fix] = .1
    for j = 1:dim_prod
        AT[i_fix,j] = 0
    end
    AT[i_fix, i_fix] = 1

    # Solve linear system for distribution
    gg = AT \ b
    g_sum = gg'*ones(dim_prod, 1) * prod(partitions)
    gg = gg ./ g_sum
    return vec(gg)
end

# Iterative method for ill-conditioned A
function solve_kfe(A::SparseMatrixCSC{Float64, Int64}, g0::Vector{Float64},
                   weight_mat::SparseMatrixCSC{Float64, Int64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) != size(A)
        error("Dimension of weight matrix is incorrect.")
    end
    dim_size = size(A, 1)
    gg = vec(g0) # stack distribution matrix into vector

    # Solve linear system
    for ikfe = 1:maxit_kfe
        gg_tilde = weight_mat * gg # weight distribution points by their measure across wealth

        # 0 = A'g -> g = A'g + g -> (I - A')g = g, but possible singularity so, given scalar Delta_KFE,
        # 0 = A' Delta_KFE g -> (I - Delta_KFE A') g = g
        gg1_tilde = (speye(Float64, dim_size) - Δ_kfe * A') \ gg_tilde
        gg1_tilde = gg1_tilde ./ sum(gg1_tilde) # normalize
        gg1 = weight_mat \ gg1_tilde # undo integration over wealth

        # Check iteration for convergence
        err_kfe = maximum(abs.(gg1 - gg))
        if err_kfe < tol_kfe
            break
        end

        gg = gg1
    end

    return vec(gg)
end
# Other iterative methods when weight_mat is (1) not sparse,
# (2) a matrix that integrates out non-wealth variables,
# (3) a vector of differences over the wealth grid
function solve_kfe(A::SparseMatrixCSC{Float64, Int64}, g0::Vector{Float64},
                   weight_mat::Matrix{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    if size(weight_mat) == size(A)
        return solve_kfe(A, g0, sparse(weight_mat); maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    else
        weight_mat = spdiagm(vec(weight_mat), 0)
        return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
    end
end
function solve_kfe(A::SparseMatrixCSC{Float64, Int64}, g0::Vector{Float64},
                   grid::Vector{Float64};
                   maxit_kfe::Int64 = 1000, tol_kfe::Float64 = 1e-12,
                   Δ_kfe::Float64 = 1e6)
    reps = size(A, 1) / length(grid)
    if Float64(reps) != Float64(Int64(reps))
        error("Dimension of A does not match dimensions of the grid.")
    end
    diagonal = zeros(size(A, 1))
    for i = 1:reps
        reps[1 + length(grid) * (i - 1): length(grid) * i] = grid
    end
    weight_mat = spdiagm(diagonal, 0)
    return solve_kfe(A, g0, weight_mat; maxit_kfe = maxit_kfe, tol_kfe = tol_kfe, Δ_kfe = Δ_kfe)
end
