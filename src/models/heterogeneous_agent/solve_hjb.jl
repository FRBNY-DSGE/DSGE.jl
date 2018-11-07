# This file holds functions that compute the upwinding scheme
# for the value function, given inputs for the savings over
# the state space, and solves the HJB equation,
# given the standard form: ρV = u + AV.

# Upwind step for value function
# Exact upwind, assumes If, Ib, and I0 do not conflict
function upwind_value_function(dVf::Matrix{T}, dVb::Matrix{T}, dV0::Matrix{T},
                               sf::Matrix{T}, sb::Matrix{T}) where {T<:Number}
    If = sf .> 0 # positive drift -> forward diff
    Ib = sb .< 0 # negative drift -> backward dift
    I0 = (1 .- If .- Ib)
    V  =  dVf .* If + dVb .* Ib + dV0 .* I0 # Compute upwinded value function difference matrix
    return V, If, Ib, I0
end

# Inexact upwind for vhen we know exact form of value function (e.g. CRRA w/ no labor disutility)
# and for when we are unsure if If, Ib, and I0 will not uniquely identify the upwind scheme.
# This function has not been verified to work yet since we have not used a model with it yet.
# This was based off the following upwind_value_function from OneAssetHANK
function upwind_value_function(dVf::Matrix{T}, dVb::Matrix{T}, dV0::Matrix{T},
                               cf::Matrix{T}, cb::Matrix{T}, c0::Matrix{T},
                               sf::Matrix{T}, sb::Matrix{T}; reverse_sign::Float64 = -1e12) where {T<:Number}
    Vf = (cf .> 0) .* (sf .* dVf) + (cf .<= 0) * (reverse_sign)
    Vb = (cb .> 0) .* (sb .* dVb) + (cb .<= 0) * (reverse_sign)
    V0 = (c0 .<= 0) * (reverse_sign)

    I_neither = (1 - (sf .> 0)) .* (1 - (sb .< 0)) # exactly zero drift
    I_unique = (sb .< 0) .* (1 - (sf .> 0)) + (1 - (sb .< 0)) .* (sf .> 0) # just one direction
    I_both = (sb .< 0) .* (sf .> 0) # both positive and negative

    Ib = I_unique .* (sb .< 0) .* (Vb .> V0) + I_both .* (Vb .== max.(max.(Vb, Vf), V0))
    If = I_unique .* (sf .> 0) .* (Vf .> V0) + I_both .* (Vf .== max.(max.(Vb, Vf), V0))
    I0 = I_neither + (1 - I_neither) .* (V0 .== max.(max.(Vb, Vf), V0))
    I0 = 1 - Ib - If
    V = dVf .* If + dVb .* Ib + dV0 .* I0

    return V, If, Ib, I0
end

# Inexact upwind, assumes positive consumption drift corresponds to positive value function drift,
# written based on model using CRRA utility with labor disutility.
# The upwind is inexact b/c it accomodates for when the savings difference indicates
# we should use both a forward and a backward difference at the same location in the upwind scheme.
# Vaf - value function without labor disutility applied to forward consumption diff
# Vf - value function w/labor disutility, etc., applied to forward consumption diff, forward hours diff, etc.
function upwind_value_function(Vaf::Matrix{T}, Vab::Matrix{T}, Va0::Matrix{T},
                               Vf_orig::Matrix{T}, Vb_orig::Matrix{T}, V0_orig::Matrix{T},
                               cf::Matrix{T}, cb::Matrix{T}, c0::Matrix{T},
                               sf::Matrix{T}, sb::Matrix{T}; reverse_sign::Float64 = -1e12) where {T<:Number}
    # Approximate forward difference in value function? Not entirely sure.
    # The Vf_orig + sf .* Vaf seems to be an approximation of a differencing of the
    # flow utility + ∂_a V * da_t
    # term in the HJB equation. See the determining of indices for another indication
    # why this seems to be some kind of approximation method.
    Vf = (cf .> 0) .* (Vf_orig + sf .* Vaf) + (cf .<= 0) * (reverse_sign)
    Vb = (cb .> 0) .* (Vb_orig + sb .* Vab) + (cb .<= 0) * (reverse_sign)
    V0 = (c0 .> 0) .* V0_orig + (c0 .<= 0) * (reverse_sign)

    I_neither = (1 - (sf .> 0)) .* (1 - (sb .< 0)) # exactly zero drift at this point
    I_unique = (sb .< 0) .* (1 - (sf .> 0)) + (1 - (sb .< 0)) .* (sf .> 0) # just one direction
    I_both = (sb .< 0) .* (sf .> 0) # both positive and negative

    # Determine indices with forward, backward, and no drift
    # Observe that Vb .> V0 is how we are determining when to consider I_unique or not
    Ib = I_unique .* (sb .< 0) .* (Vb .> V0) + I_both .* (Vb .== max.(max.(Vb, Vf), V0))
    If = I_unique .* (sf .> 0) .* (Vf .> V0) + I_both .* (Vf .== max.(max.(Vb, Vf), V0))
    I0 = I_neither + (1 - I_neither) .* (V0 .== max.(max.(Vb, Vf), V0))
    I0 = 1 - Ib - If

    # Compute upwinded value function difference (in wealth direction) matrix
    V = Vaf .* If + Vab .* Ib + Va0 .* I0

    return V, If, Ib, I0
end

# Constructs A matrix that applies to the HJB
# A_switch - summarizes all information due to non-wealth state variables, e.g. income shocks
function upwind_matrix(A_switch::SparseMatrixCSC{R, Int64},
                       sf::Matrix{S}, sb::Matrix{S},
                       f_diff::Matrix{T}, b_diff::Matrix{T},
                       wealth_dim::Int, other_dims::Int;
                       exact::Bool = true,
                       If::Matrix{Int64} = Matrix{Int64}(0, 0),
                       Ib::Matrix{Int64} = Matrix{Int64}(0, 0),
                       I0::Matrix{Int64} = Matrix{Int64}(0, 0)) where {R<:Number, S<:Number, T<:Number}
    if size(sf, 1) != wealth_dim
        error("Dimension of wealth grid, wealth_dim, is incorrect.")
    elseif size(f_diff) != size(sf) || size(b_diff) != size(sf)
        error("Dimension of diff, differenced state space grid, along the wealth dimension is wrong.")
    end

    # Compute diagonals of first-order differenced V matrix
    if exact
        # use min/max b/c assumes we set I0 entries to 1 only if sf/sb exactly zero
        X = -min.(sb, 0) ./ b_diff # If savings drift < 0, re-weight by grid partition
                                   # (to approx a derivative), and set to negative to get the
                                   # right direction b/c lower diagonal (doing a backward
                                   # difference here so set positive)
        Y = -max.(sf, 0) ./ f_diff + min.(sb, 0) ./ b_diff # Same idea here
        Z = max.(sf, 0) ./ f_diff # Same idea here
    else
        X = -Ib .* sb ./ b_diff
        Y = -If .* sf ./ f_diff + Ib .* sb ./ b_diff
        Z =  If .* sf ./ f_diff
    end

    total_dims = wealth_dim * other_dims
    centdiag = vec(Y) # Y is our main diagonal
    updiag = zeros(R, total_dims - 1)
    lowdiag = zeros(R, total_dims - 1)
    for j = 1:other_dims
        for i = 1:wealth_dim
            # When i == wealth_dim, leave entry as zero b/c matrix should be block
            # [nonzero zeros; zeros nonzero] in the case of a wealth_dim x 2 state space
            # Also, this way, we don't save Z[end, end]
            # or X[1, 1] since X and Z are lower/upper diagonals
            if i < wealth_dim
                updiag[wealth_dim * (j-1) + i] = Z[i, j]
                lowdiag[wealth_dim * (j-1) + i] = X[i + 1, j]
            end
        end
    end

    AA = spdiagm((lowdiag, centdiag, updiag), (-1, 0, 1)) # AA * v ≈ ∂_a v * (time derivative of wealth)
    A = AA + A_switch # main A matrix in standard representation of HACT steady state
    return A
end

# Same function as above but with a constant difference in the state space grid
function upwind_matrix(A_switch::SparseMatrixCSC{R, Int64},
                       sf::Matrix{S}, sb::Matrix{S},
                       diff::T, wealth_dim::Int, other_dims::Int;
                       exact::Bool = true,
                       If::Matrix{Int64} = Matrix{Int64}(undef, 0, 0),
                       Ib::Matrix{Int64} = Matrix{Int64}(undef, 0, 0),
                       I0::Matrix{Int64} = Matrix{Int64}(undef, 0, 0)) where {R<:Number, S<:Number, T<:Number}
    if size(sf, 1) != wealth_dim
        error("Dimension of wealth grid is incorrect.")
    end

    # Compute diagonals of first-order differenced V matrix
    if exact
        # use min/max b/c assumes we set I0 entries to 1 only if sf/sb exactly zero
        X = -min.(sb, 0) / diff
        Y = -max.(sf, 0) / diff + min.(sb, 0) / diff
        Z =  max.(sf, 0) / diff
    else
        # Need to multiply with Ib, If to be consistent with how we set the entries of I0
        X = -Ib .* sb ./ diff
        Y = -If .* sf ./ diff + Ib .* sb ./ diff
        Z =  If .* sf ./ diff
    end
    total_dims = wealth_dim * other_dims
    centdiag = vec(Y)
    updiag = zeros(R, total_dims - 1)
    lowdiag = zeros(R, total_dims - 1)
    for j = 1:other_dims
        for i = 1:wealth_dim
            if i < wealth_dim
                updiag[wealth_dim * (j-1) + i] = Z[i, j]
                lowdiag[wealth_dim * (j-1) + i] = X[i + 1, j]
            end
        end
    end

    AA = SparseArrays.spdiagm(0 => centdiag, 1 => updiag, -1 => lowdiag)
    A = AA + A_switch
    return A
end

# Solves hjb via left divide
# Δ_HJB is used to make the left-divide more numerically stable. Derivation is as follows:
# (ρ - A) * V = u
# (ρ - A) * V + V/Δ_HJB = u + V/Δ_HJB
# Since Δ_HJB is large, V/Δ_HJB ≈ 0 but also adds a small perturbation to reflect
# that we are approximating V, so our vector u may not be exactly correct. Re-arranging yields
# the form ((1/Δ_HJB + ρ) - A) V = u + V/Δ_HJB
function solve_hjb(A::SparseMatrixCSC{R, Int64}, ρ::S,
                   Δ_HJB::S, u::Matrix{T}, V::Matrix{T}) where {R<:Number, S<:Number, T<:Number}
    # Stack u and V
    return solve_hjb(A, ρ, Δ_HJB, vec(u), vec(V))
end

function solve_hjb(A::SparseMatrixCSC{R, Int64}, ρ::S,
                   Δ_HJB::S, u::Vector{T}, V::Vector{T}) where {R<:Number, S<:Number, T<:Number}
    B = (1/Δ_HJB + ρ) * SparseMatrixCSC(I, size(A, 1), size(A, 1)) - A
    b = u + V / Δ_HJB
    return vec(B \ b)
end

# Wrapper function combining upwind_value_function and upwind_matrix
# since these do not require intermediate steps. However, b/c
# solve_hjb needs flow utility, which is model-specific,
# we do not wrap all of these steps into one function
function upwind(A_switch::SparseMatrixCSC{R, Int64},
                Vaf::Matrix{S}, Vab::Matrix{S}, Va0::Matrix{S},
                Vf::Matrix{S}, Vb::Matrix{S}, V0::Matrix{S},
                cf::Matrix{S}, cb::Matrix{S}, c0::Matrix{S},
                sf::Matrix{S}, sb::Matrix{S},
                f_diff::Matrix{T}, b_diff::Matrix{T},
                wealth_dim::Int, other_dims::Int;
                reverse_sign = -1e12, exact::Bool = true) where {R<:Number, S<:Number, T<:Number}
    if exact
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, sf, sb)
        A = upwind_matrix(A_switch, sf, sb, f_diff, b_diff, wealth_dim, other_dims;
                      exact = exact)
    else
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, Vf, Vb, V0, cf, cb, c0,
                                                  sf, sb; reverse_sign = reverse_sign)
        A = upwind_matrix(A_switch, sf, sb, f_diff, b_diff, wealth_dim, other_dims;
                      exact = exact, If = If, Ib = Ib, I0 = I0)
    end
    return dV_upwind, If, Ib, I0, A
end

function upwind(A_switch::SparseMatrixCSC{R, Int64},
                Vaf::Matrix{S}, Vab::Matrix{S}, Va0::Matrix{S},
                sf::Matrix{S}, sb::Matrix{S},
                f_diff::Matrix{T}, b_diff::Matrix{T},
                wealth_dim::Int, other_dims::Int,
                cf::Matrix{S} = Matrix{S}(0, 0),
                cb::Matrix{S} = Matrix{S}(0, 0),
                c0::Matrix{S} = Matrix{S}(0, 0);
                exact::Bool = true,
                reverse_sign::Float64 = -1e12) where {R<:Number, S<:Number, T<:Number}

    if exact
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, sf, sb)
        A = upwind_matrix(A_switch, sf, sb, f_diff, b_diff, wealth_dim, other_dims;
                      exact = exact)
    else
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, cf, cb, c0, sf, sb;
                                                      reverse_sign = reverse_sign)
        A = upwind_matrix(A_switch, sf, sb, f_diff, b_diff, wealth_dim, other_dims;
                      exact = exact, If = If, Ib = Ib, I0 = I0)
    end

    return dV_upwind, If, Ib, I0, A
end

function upwind(A_switch::SparseMatrixCSC{R, Int64},
                Vaf::Matrix{S}, Vab::Matrix{S}, Va0::Matrix{S},
                Vf::Matrix{S}, Vb::Matrix{S}, V0::Matrix{S},
                cf::Matrix{S}, cb::Matrix{S}, c0::Matrix{S},
                sf::Matrix{S}, sb::Matrix{S},
                diff::T, wealth_dim::Int, other_dims::Int;
                reverse_sign = -1e12, exact::Bool = trie) where {R<:Number, S<:Number, T<:Number}
    if exact
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, sf, sb)
        A = upwind_matrix(A_switch, sf, sb, diff, wealth_dim, other_dims;
                      exact = exact)
    else
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, Vf, Vb, V0, cf, cb, c0,
                                                  sf, sb; reverse_sign = reverse_sign)
        A = upwind_matrix(A_switch, sf, sb, diff, wealth_dim, other_dims;
                          exact = exact, If = If, Ib = Ib, I0 = I0)
    end
    return dV_upwind, If, Ib, I0, A
end

function upwind(A_switch::SparseMatrixCSC{R, Int64},
                Vaf::Matrix{S}, Vab::Matrix{S}, Va0::Matrix{S},
                sf::Matrix{S}, sb::Matrix{S},
                diff::T, wealth_dim::Int, other_dims::Int;
                exact::Bool = true,
                reverse_sign::Float64 = -1e12) where {R<:Number, S<:Number, T<:Number}
    if exact
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, sf, sb)
        A = upwind_matrix(A_switch, sf, sb, diff, wealth_dim, other_dims;
                      exact = exact)
    else
        dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, cf, cb, c0, sf, sb;
                                                      reverse_sign = reverse_sign)
        A = upwind_matrix(A_switch, sf, sb, diff, wealth_dim, other_dims;
                      exact = exact, If = If, Ib = Ib, I0 = I0)
    end

    return dV_upwind, If, Ib, I0, A
end