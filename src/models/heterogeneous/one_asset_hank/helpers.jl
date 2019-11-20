# This file contains additional helper functions for computing the steady state.

# Discretize a state space dimension of I gridpoints
# Where agridparam is the "bending" coefficient of the grid
# i.e. grid_param = 1 implies a uniform grid
# And grid_min/grid_max are the start and endpoints of the grid
@inline function construct_asset_grid(I::Int64, grid_param::Int64, grid_min::Float64, grid_max::Float64)
    a  = collect(range(0, stop=1, length=I))
    a  = a .^ (1 / grid_param)
    a  = grid_min .+ (grid_max - grid_min) * a
    return a
end

# Provide the income and asset grids
@inline function initialize_diff_grids(a::Vector{Float64}, I::Int64, J::Int64)
    daf    = similar(a) # forward difference for a
    dab    = similar(a) # backward difference for a
    adelta = similar(a) # size of differences

    for i=1:I
        # Create a grid of lengths of overlapping intervals in a dimension.
        # The purpose is generally to compute Riemann integrals
        # by taking midpoint Riemann sums and dividing by two to adjust for
        # or average out the added Lebesgue measure given by using overlapping intervals.
        daf[i] = (i==I) ? a[end] - a[end-1] : a[i+1] - a[i]
        dab[i] = (i==1) ? a[2] - a[1]       : a[i] - a[i-1]

        # Ex:    adelta[2] =    (a_3 - a_1) / 2; if we integrate f, then
        # f[2] * adelta[2] = f_2(a_3 - a_1) / 2
        if i==1
            adelta[1]   = 0.5 * daf[1]
        elseif i==I
            adelta[end] = 0.5 * daf[end-1]
        else
            adelta[i]   = 0.5 * (daf[i-1] + daf[i])
        end
    end
    azdelta = repeat(adelta, J)

    return daf, dab, azdelta
end

# P is the Markov transition matrix
# n_income_states is the number of income states (the number of discrete states in the distribution)
# iter_num is the number of iterations willing to be accepted for convergence
@inline function compute_stationary_income_distribution(P::Matrix{Float64}, n_income_states::Int64;
                                                iter_num::Int64 = 50)
    Pt = P'
    g_z = fill(1/n_income_states, n_income_states)
    for n = 1:iter_num
        g_z_new = (speye(n_income_states) - Pt * 1000)\g_z
        diff    = maximum(abs.(g_z_new - g_z))
        if diff < 1e-5
            break
        end
        g_z = g_z_new
    end
    return g_z
end

# initial_ygrid is the unscaled income values for each income state
# income_distr is the stationary income distribution
# meanlabeff is the mean labor efficiency value that scales z such that
# The expected value of income with respect to income_distr is meanlabeff
# n_gridpoints should be the same as the number of grid points in the asset process
@inline function construct_labor_income_grid(initial_ygrid::Vector{Float64},
                                     income_distr::Vector{Float64},
                                     meanlabeff::Float64, n_gridpoints::Int64)
    z       = exp.(initial_ygrid)
    z_bar   = dot(z, income_distr)
    z  	    = (meanlabeff/z_bar) .* z
    zz      = ones(n_gridpoints, 1) * z'
    return zz
end

@inline function construct_household_problem_functions(V::Matrix{S}, w::T,
                                                       coefrra::R, frisch::R,
                                                       labtax::R,
                                                       labdisutil::R) where {R<:AbstractFloat,T<:Real,S<:Number}

    @inline function util(c::U, h::U) where {U<:Number}
        f(x::U) = coefrra == 1.0 ? log(x) : x^(1-coefrra) / (1-coefrra)
        return f(c) - labdisutil * (h ^ (1 + 1/frisch)/(1 + 1/frisch))
    end

    @inline income(h::U, z::Float64, profshare::V, lumptransfer::V, r::T,
                   a::Float64) where {T<:Number,U<:Number,V<:Number} = h * z * w * (1 - labtax) + lumptransfer + profshare + r * a

    @inline labor(z::U, val::V) where {U<:Number,V<:Number} = (z * w * (1 - labtax) * val / labdisutil) ^ frisch

    return util, income, labor
end

@inline function construct_initial_diff_matrices(V::Matrix{T},
                                                 Vaf::Matrix{T}, Vab::Matrix{T},
                                                 income::Function, labor::Function,
                                                 h::Matrix{U}, h0::Matrix{U},
                                                 zz::Matrix{S}, profshare::Matrix{T},
                                                 lumptransfer::T,
                                                 amax::S, amin::S, coefrra::S, r::S,
                                                 daf::Vector{S}, dab::Vector{S},
                                                 maxhours::S) where {S<:Number,T<:Number,U<:Number}
    I,J = size(V)
    cf  = similar(V)
    hf  = similar(V)
    cb  = similar(V)
    hb  = similar(V)

    for j=1:J, i=1:I
        if i==I
            Vaf[end, j] = income(h[end, j], zz[end, j], profshare[end, j], lumptransfer, r, amax) ^ (-coefrra)
            Vab[i, j]   = (V[i, j] - V[i-1, j]) / dab[i]
        elseif i==1
            Vaf[i, j]   = (V[i+1, j] - V[i, j]) / daf[i]
            Vab[1, j]   = income(h0[1, j], zz[1, j], profshare[1, j], lumptransfer, r, amin) ^ (-coefrra)
        else
            Vaf[i, j]   = (V[i+1, j] - V[i, j]) / daf[i]
            Vab[i, j]   = (V[i, j] - V[i-1, j]) / dab[i]
        end
        cf[i,j] = Vaf[i,j] ^ (-1 / coefrra)
        cb[i,j] = Vab[i,j] ^ (-1 / coefrra)

        hf[i,j] = min(norm(labor(zz[i,j], Vaf[i,j])), maxhours)
        hb[i,j] = min(norm(labor(zz[i,j], Vab[i,j])), maxhours)
    end

    return Vaf, Vab, cf, hf, cb, hb
end

# For initialization
@inline function calculate_ss_equil_vars(zz::Matrix{Float64}, m_ss::Float64, meanlabeff::Float64,
                                 lumptransferpc::Float64, govbondtarget::Float64)

    N_ss         = complex(1/3) # steady state hours: so that quarterly GDP = 1 in s.s
    Y_ss         = complex(1.)
    B_ss         = govbondtarget * Y_ss
    profit_ss    = complex((1 - m_ss) * Y_ss)
    profshare    = zz / meanlabeff * profit_ss
    lumptransfer = complex(lumptransferpc * Y_ss)

    return N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer
end

@inline function calculate_ss_equil_vars(zz::Matrix{Float64},
                                 h::Matrix{ComplexF64},
                                 g::Vector{ComplexF64}, azdelta::Vector{Float64},
                                 aa::Matrix{Float64}, m_ss::Float64,
                                 meanlabeff::Float64, lumptransferpc::Float64,
                                 govbondtarget::Float64)

    # equilibrium objects
    Y_ss = N_ss  = sum(vec(zz) .* vec(h) .* g .* azdelta)
    Y_ss         = N_ss
    B_ss         = sum(g .* vec(aa) .* azdelta)
    profit_ss    = (1 - m_ss) * Y_ss
    profshare    = zz / meanlabeff * profit_ss
    lumptransfer = lumptransferpc * Y_ss
    bond_err     = B_ss / Y_ss - govbondtarget

    return N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer, bond_err
end

@inline function hours_iteration(income::Function, labor::Function,
                                 zz::Matrix{Float64},
                                 profshare::Matrix{T},
                                 lumptransfer::T,
                                 aa::Matrix{Float64},
                                 coefrra::Float64, r::S,
                                 cf::Matrix{T}, hf::Matrix{T},
                                 cb::Matrix{T}, hb::Matrix{T},
                                 c0::Matrix{T}, h0::Matrix{U},
                                 maxhours::Float64,
                                 niter_hours::Int64) where {S<:Number,T<:Number,U<:Number}
    I, J = size(zz)
    for ih = 1:niter_hours
        for j=1:J, i=1:I
            if i==I
                cf[end, j] = income(hf[end, j], zz[end, j], profshare[end, j], lumptransfer, r, aa[end, j])
                hf[end, j] = labor(zz[end, j], cf[end, j] ^ (-coefrra))
                hf[end, j] = min(norm(hf[end, j]), maxhours)
            elseif i==1
                cb[1, j] = income(hb[1, j], zz[1, j], profshare[1, j], lumptransfer, r, aa[1, j])
                hb[1, j] = labor(zz[1, j], cb[1, j] ^ (-coefrra))
                hb[1, j] = min(norm(hb[1, j]), maxhours)
            end
            c0[i, j] = income(h0[i, j], zz[i, j], profshare[i, j], lumptransfer, r, aa[i, j])
            h0[i, j] = labor(zz[i, j], c0[i, j]^(-coefrra))
            h0[i, j] = min(norm(h0[i, j]), maxhours)
        end
    end
    return cf, hf, cb, hb, c0, h0
end

# For compute_steady_state; choose upwinding direction
@inline function upwind(ρ::Float64, V::Matrix{T}, args...;
                        Δ_HJB::Float64 = 1e6) where {T<:Number}
    A, u, h, c, s = upwind(args...)
    I, J = size(u)
    B    = (1 / Δ_HJB + ρ) * speye(T, I*J) - A
    b    = reshape(u, I*J) + reshape(V, I*J) / Δ_HJB
    V    = reshape(B \ b, I, J)
    return V, A, u, h, c, s
end

# For equilibrium_conditions; choose upwinding direction.
@inline function upwind(util::Function,
                        A_switch::SparseMatrixCSC{S, Int64},
                        cf::Matrix{T}, cb::Matrix{T},
                        c0::Matrix{T}, hf::Matrix{T},
                        hb::Matrix{T}, h0::Matrix{T},
                        sf::Matrix{T}, sb::Matrix{T},
                        Vaf::Matrix{T}, Vab::Matrix{T},
                        daf::Vector{Float64},
                        dab::Vector{Float64}) where {T<:Number, S<:Number}

    I,J = size(sb)
    h   = similar(sb)
    c   = similar(sb)
    s   = similar(sb)
    u   = similar(sb)
    X   = similar(sb)
    Z   = similar(sb)
    Y   = similar(sb)

    for i in eachindex(sb)

        Vf = (cf[i] > 0) * (util(cf[i], hf[i]) + sf[i] * Vaf[i]) + (cf[i] <= 0) * (-1e12)
        Vb = (cb[i] > 0) * (util(cb[i], hb[i]) + sb[i] * Vab[i]) + (cb[i] <= 0) * (-1e12)
        V0 = (c0[i] > 0) * util(c0[i], h0[i]) + (c0[i] <= 0) * (-1e12)

        Iunique = (sb[i] < 0) * (1 - (sf[i] > 0)) + (1 - (sb[i] < 0)) * (sf[i] > 0)
        Iboth = (sb[i] < 0) * (sf[i] > 0)

        Ib = Iunique * (sb[i] < 0) * (Vb > V0) + Iboth * (Vb == max(max(Vb, Vf), V0))
        If = Iunique * (sf[i] > 0) * (Vf > V0) + Iboth * (Vf == max(max(Vb, Vf), V0))
        I0 = 1 - Ib - If

        h[i] = hf[i] * If + hb[i] * Ib + h0[i] * I0
        c[i] = cf[i] * If + cb[i] * Ib + c0[i] * I0
        s[i] = sf[i] * If + sb[i] * Ib
        u[i] = util(c[i], h[i])

        # Construct A matrix
        X[i] = -Ib * sb[i] / dab[((i-1) % I) + 1]
        Z[i] =  If * sf[i] / daf[((i-1) % I) + 1]
        Y[i] = -Z[i] - X[i]
    end

    # R: Pretty sure the indexing for this ought be the same as in the KrusellSmith model
    X[1,:] .= T == ComplexF64 ? complex(0.) : 0.
    Z[I,:] .= T == ComplexF64 ? complex(0.) : 0.

    A = spdiagm(-1 => reshape(X,I*J)[2:I*J], 0 => reshape(Y,I*J), 1 => reshape(Z,I*J)[1:I*J-1]) + A_switch

    return A, u, h, c, s
end


# Using the market clearing condition on bonds to determine whether or not
# an equilibrium has been reached
@inline function check_bond_market_clearing(bond_err::ComplexF64, crit_S::Float64,
                                            r::Float64, r_min::Float64, r_max::Float64,
                                            r_ρ::Float64, ρ_min::Float64, ρ_max::Float64,
                                            iter_r::Bool, iter_ρ::Bool)
    clearing_condition = false
    # Using the market clearing condition on bonds to determine whether or not
    # an equilibrium has been reached
    if abs(bond_err) > crit_S
        if bond_err > 0
            if iter_r
                r_max  = r
                r      = 0.5 * (r + r_min)
            elseif iter_ρ
                ρ_min = r_ρ
                r_ρ   = 0.5 * (r_ρ + ρ_max)
            end
        else
            if iter_r
                r_min  = r
                r      = 0.5 * (r + r_max)
            elseif iter_ρ
                ρ_max = r_ρ
                r_ρ   = 0.5 * (r_ρ + ρ_min)
            end
        end
    else
        clearing_condition = true
    end
    return r, r_min, r_max, r_ρ, ρ_min, ρ_max, clearing_condition
end