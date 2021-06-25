"""
```
@inline function solve_hjb(V::Union{Array{R},Array{S}}, a_lb::R, γ::R, ddeath::R, pam::R,
                           trans::R, ξ::R, τ_I::R, aggZ::Union{R,S}, w::Union{R,S},
                           r_b_vec::Union{Vector{R},Vector{S}}, y::Vector{U}, a::Vector{R},
                           b::Vector{R}, cost, util, deposit; permanent::Bool=false,
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, U<:Number}
```

Solves Hamilton-Jacobi-Bellman equation.
"""
@inline function solve_hjb(V::Union{Array{R},Array{S}}, a_lb::R, γ::R, ddeath::R, pam::R,
                           trans::R, ξ::R, τ_I::R, aggZ::Union{R,S}, w::Union{R,S},
                           r_b_vec::Union{Vector{R},Vector{S}}, y::Vector{U}, a::Vector{R},
                           b::Vector{R}, cost, util, deposit; permanent::Bool=false,
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, U<:Number}
    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------
    I, J, N = size(V)
    c, s, d = similar(V), similar(V), similar(V)

    perm_c = permanent ? ddeath * pam - aggZ : ddeath * pam
    c0_c   = ((1-ξ) - τ_I) * w

    for i=1:I, j=1:J, n=1:N

        c0 = c0_c * real(y[n]) + b[i] * (r_b_vec[i] + perm_c) + trans

        # ---- Liquid Assets, Forward + Backward Difference ----
        VaF = (j==J) ? 0.0 : max((V[i,j+1,n] - V[i,j,n]) / (a[j+1]-a[j]), Vamin)
        VaB = (j==1) ? 0.0 : max((V[i,j,n] - V[i,j-1,n]) / (a[j]-a[j-1]), Vamin)

        # ---- Illiquid Assets, Forward + Backward Difference ----
        VbF = (i==I) ? 0.0 : max((V[i+1,j,n] - V[i,j,n]) / (b[i+1]-b[i]), Vbmin)
        VbB = (i==1) ? 0.0 : max((V[i,j,n] - V[i-1,j,n]) / (b[i]-b[i-1]), Vbmin)

        # Decisions conditional on a particular direction of derivative
        cF = (i==I) ? 0.0 : VbF ^ (-1 / γ)
        cB = (i==1) ? c0  : VbB ^ (-1 / γ)

        sF = (i==I) ? 0.0 : c0 - cF
        sB = (i==1) ? 0.0 : c0 - cB

        if (i==1 && j > 1) VbB = util(cB) end

        #----------------------------------------------------------------
        # Consumption & Savings Decision
        #----------------------------------------------------------------
        Hc0 =                  util(c0)
        HcF = (i==I) ? -1e12 : util(cF) + VbF * sF
        HcB =                  util(cB) + VbB * sB

        validF = (sF > 0)
        validB = (sB < 0)

        IcF = validF * max(!validB, (HcF >= HcB)) * (HcF >= Hc0)
        IcB = validB * max(!validF, (HcB >= HcF)) * (HcB >= Hc0)
        Ic0 = 1 - IcF - IcB

        c[i,j,n] = IcF * cF + IcB * cB + Ic0 * c0
        s[i,j,n] = IcF * sF + IcB * sB

        #----------------------------------------------------------------
        # Deposit Decision
        #----------------------------------------------------------------
        dFB = (i == 1 || j == J) ? 0.0 : deposit(VaF, VbB, max(a[j], a_lb))
        dBF = (i == I || j == 1) ? 0.0 : deposit(VaB, VbF, max(a[j], a_lb))
        dBB = (j == 1)           ? 0.0 : deposit(VaB, VbB, max(a[j], a_lb))

        HdFB = (i == 1 || j == J) ? -1e12 : VaF * dFB - VbB * (dFB + cost(dFB, a[j]))
        HdBB = (j == 1)           ? -1e12 : VaB * dBB - VbB * (dBB + cost(dBB, a[j]))
        HdBF = (i == I || j == 1) ? -1e12 : VaB * dBF - VbF * (dBF + cost(dBF, a[j]))

        validFB = (HdFB > 0) * (dFB >  0)
        validBB = (HdBB > 0) * (dBB <= 0) * (dBB >  -cost(dBB, a[j]))
        validBF = (HdBF > 0) *              (dBF <= -cost(dBF, a[j]))

        IcFB = max(!validBF, HdFB >= HdBF) * validFB * max(!validBB, HdFB >= HdBB)
        IcBF = max(!validFB, HdBF >= HdFB) * validBF * max(!validBB, HdBF >= HdBB)
        IcBB = max(!validFB, HdBB >= HdFB) * validBB * max(!validBF, HdBB >= HdBF)

        d[i,j,n] = IcFB * dFB + IcBF * dBF + IcBB * dBB
    end
    return c, s, d
end

"""
```
@inline function solve_hjb_lite(V::Union{Array{R},Array{S}}, a_lb::R, γ::R,
                           permanent::Bool, ddeath::R, pam::R, aggZ::Union{R,S}, ξ::R, τ_I::R,
                           w::Union{R,S}, trans::R, r_b::Union{R,S}, r_b_borr::Union{R,S},
                           y::U, a::Vector{R}, b::Vector{R}, cost, util, deposit;
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, U<:Number}
```
Solves Hamilton-Jacobi-Bellman equation.
"""
@inline function solve_hjb_lite(V::Union{Array{R},Array{S}}, a_lb::R, γ::R,
                           permanent::Bool, ddeath::R, pam::R, aggZ::Union{R,S}, ξ::R, τ_I::R,
                           w::Union{R,S}, trans::R, r_b::Union{R,S}, r_b_borr::Union{R,S},
                           y::U, a::Vector{R}, b::Vector{R}, cost, util, deposit;
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, U<:Number}
    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------
    I, J = size(V)

    perm_c = permanent ? ddeath * pam - aggZ : ddeath * pam
    c0_c   = ((1-ξ) - τ_I) * w

    c, s, d = similar(V), similar(V), similar(V)
    for i=1:I, j=1:J
        c0 = c0_c * real(y) + b[i] * (((b[i] >= 0) ? r_b : r_b_borr) + perm_c) + trans

        # ---- Liquid Assets, Forward + Backward Difference ----
        VaF = (j==J) ? 0.0 : max((V[i,j+1] - V[i,j]) / (a[j+1]-a[j]), Vamin)
        VaB = (j==1) ? 0.0 : max((V[i,j] - V[i,j-1]) / (a[j]-a[j-1]), Vamin)

        # ---- Illiquid Assets, Forward + Backward Difference ----
        VbF = (i==I) ? 0.0 : max((V[i+1,j] - V[i,j]) / (b[i+1]-b[i]), Vbmin)
        VbB = (i==1) ? 0.0 : max((V[i,j] - V[i-1,j]) / (b[i]-b[i-1]), Vbmin)

        # Decisions conditional on a particular direction of derivative
        cF = (i==I) ? 0.0 : VbF ^ (-1 / γ)
        cB = (i==1) ? c0  : VbB ^ (-1 / γ)

        sF = (i==I) ? 0.0 : c0 - cF
        sB = (i==1) ? 0.0 : c0 - cB

        if (i==1 && j > 1) VbB = util(cB) end

        #----------------------------------------------------------------
        # Consumption & Savings Decision
        #----------------------------------------------------------------
        Hc0 =                  util(c0)
        HcF = (i==I) ? -1e12 : util(cF) + VbF * sF
        HcB =                  util(cB) + VbB * sB

        validF = (sF > 0)
        validB = (sB < 0)

        IcF = validF * max(!validB, (HcF >= HcB)) * (HcF >= Hc0)
        IcB = validB * max(!validF, (HcB >= HcF)) * (HcB >= Hc0)
        Ic0 = 1 - IcF - IcB

        c[i,j] = IcF * cF + IcB * cB + Ic0 * c0
        s[i,j] = IcF * sF + IcB * sB

        #----------------------------------------------------------------
        # Deposit Decision
        #----------------------------------------------------------------
        dFB = (i == 1 || j == J) ? 0.0 : deposit(VaF, VbB, max(a[j], a_lb))
        dBF = (i == I || j == 1) ? 0.0 : deposit(VaB, VbF, max(a[j], a_lb))
        dBB = (j == 1)           ? 0.0 : deposit(VaB, VbB, max(a[j], a_lb))

        HdFB = (i == 1 || j == J) ? -1e12 : VaF * dFB - VbB * (dFB + cost(dFB, a[j]))
        HdBB = (j == 1)           ? -1e12 : VaB * dBB - VbB * (dBB + cost(dBB, a[j]))
        HdBF = (i == I || j == 1) ? -1e12 : VaB * dBF - VbF * (dBF + cost(dBF, a[j]))

        validFB = (HdFB > 0) * (dFB >  0)
        validBB = (HdBB > 0) * (dBB <= 0) * (dBB >  -cost(dBB, a[j]))
        validBF = (HdBF > 0) *              (dBF <= -cost(dBF, a[j]))

        IcFB = max(!validBF, HdFB >= HdBF) * validFB * max(!validBB, HdFB >= HdBB)
        IcBF = max(!validFB, HdBF >= HdFB) * validBF * max(!validBB, HdBF >= HdBB)
        IcBB = max(!validFB, HdBB >= HdFB) * validBB * max(!validBF, HdBB >= HdBF)

        d[i,j] = IcFB * dFB + IcBF * dBF + dBB * IcBB
    end
    return c, s, d
end

@inline function solve_kfe(dab_g_tilde, b, g_inds, gg, ddeath, lambda, AT, indx, V_SS)
    my_inds = findall(x -> x != indx, 1:size(lambda,1))

    loc = findall(b .== 0)
    dab_g_tilde_inv       = vec(1.0 ./ dab_g_tilde)
    dab_g_small           = dab_g_tilde ./ dab_g_tilde[loc] * ddeath
    dab_g_small[loc]     .= 0.0
    death_process         = -ddeath * my_speye(length(gg))
    death_process[loc,:]  = vec(dab_g_small)
    gIntermediate = dab_g_tilde_inv .* ((AT * (dab_g_tilde .* gg)) +
                        sum([lambda[kk,indx] .* dab_g_tilde .* g_inds[:,kk] for kk in my_inds])) +
                                             death_process * gg
    cc = sum([lambda[indx,kk]*(V_SS[:,kk] - V_SS[:,indx]) for kk in my_inds])
    return gIntermediate, cc
end

@inline function add_diag(A::SparseMatrixCSC, c::AbstractFloat)
    for i=1:size(A,1)
        A[i,i] += c
    end
    return A
end

@inline function mul_diag(A::SparseMatrixCSC, c::AbstractFloat)
    for i=1:size(A,1)
        A[i,i] *= c
    end
    return A
end

"""
```
@inline function update_value_fn(A::SparseMatrixCSC, Vn::Union{Array{T,3},Array{R,3}},
                                 lambda::Matrix{R},   u::Union{Array{T,3},Array{R,3}},
                                 Delta::R, rrho::R, ddeath::R) where {R<:Float64, T<:Complex}
```
Update value function.
"""
@inline function update_value_fn(A::SparseMatrixCSC, Vn::Union{Array{T,3},Array{R,3}},
                                 lambda::Matrix{R},   u::Union{Array{T,3},Array{R,3}},
                                 Delta::R, rrho::R, ddeath::R) where {R<:Float64, T<:Complex}
    I, J, N = size(Vn)
    Vn_new = Array{Float64,3}(undef, I, J, N)
    for kk = 1:N
        indx_k = 1:N .!= kk
        Ak = A[1+(kk-1)*(I*J):kk*(I*J), 1+(kk-1)*(I*J):kk*(I*J)]
        Bk = add_diag(-Delta*Ak, 1 + Delta*(rrho + ddeath) - Delta*lambda[kk,kk])
        uk_stacked   = vec(u[:,:,kk])
        Vk_stacked   = vec(Vn[:,:,kk])
        Vkp_stacked  = sum(broadcast(*,  lambda[kk, indx_k]',
                                    reshape(Vn[:,:,indx_k], I*J, N-1)), dims=2)
        qk           = Delta .* uk_stacked + Vk_stacked + Delta .* Vkp_stacked
        Vn1k_stacked = Bk \ qk
        Vn_new[:,:,kk] = vec(Vn1k_stacked)
    end
    return Vn_new
end

@inline function set_grids(a, b, a_g, b_g, N)
    I   = length(b)
    J   = length(a)
    I_g = length(b_g)
    J_g = length(a_g)

    #b_grid   = permutedims(repeat(b, 1, J, N), [1 2 3])
    #a_grid   = permutedims(repeat(a, 1, I, N), [2 1 3])

    b_g_grid = permutedims(repeat(b_g, 1, J_g, N),  [1 2 3])
    a_g_grid = permutedims(repeat(a_g, 1, I_g, N),  [2 1 3])

    dbf_g_grid = Array{Float64}(undef, I_g, J_g, N)
    dbf_g_grid[1:I_g-1,:,:] = b_g_grid[2:I_g,:,:] - b_g_grid[1:I_g-1,:,:]
    dbf_g_grid[I_g,:,:] = dbf_g_grid[I_g-1,:,:]

    dbb_g_grid = Array{Float64}(undef, I_g,J_g,N)
    dbb_g_grid[2:I_g,:,:] = b_g_grid[2:I_g,:,:] - b_g_grid[1:I_g-1,:,:]
    dbb_g_grid[1,:,:] = dbb_g_grid[2,:,:]

    daf_g_grid = Array{Float64}(undef, I_g, J_g, N)
    daf_g_grid[:,1:J_g-1,:] = a_g_grid[:,2:J_g,:] - a_g_grid[:,1:J_g-1,:]
    daf_g_grid[:,J_g,:] = daf_g_grid[:,J_g-1,:]

    dab_g_grid = Array{Float64}(undef, I_g, J_g, N)
    dab_g_grid[:,2:J_g,:] = a_g_grid[:,2:J_g,:] - a_g_grid[:,1:J_g-1,:]
    dab_g_grid[:,1,:] = dab_g_grid[:,2,:]

    db_g_tilde       = 0.5*(dbb_g_grid[:,1,1] + dbf_g_grid[:,1,1])
    db_g_tilde[1]    = 0.5*dbf_g_grid[1,1,1]
    db_g_tilde[end]  = 0.5*dbb_g_grid[end,1,1]
    da_g_tilde       = 0.5*(dab_g_grid[1,:,1] + daf_g_grid[1,:,1])
    da_g_tilde[1]    = 0.5*daf_g_grid[1,1,1]
    da_g_tilde[end]  = 0.5*dab_g_grid[1,end,1]

    dab_g_tilde      = kron(da_g_tilde, db_g_tilde)
    dab_g_tilde_grid = reshape(repeat(dab_g_tilde, N, 1), I_g, J_g, N)

    return a_g_grid, b_g_grid, dab_g_tilde_grid, dab_g_tilde
end

"""
```
@inline function backward_difference(a_g::Vector{T}, b_g::Vector{T}) where {T<:Real}
```
Instantiates necessary difference vector.
"""
@inline function backward_difference(a_g::Vector{T}, b_g::Vector{T}) where {T<:Real}

    I_g, J_g = length(b_g), length(a_g)

    dbf_g = similar(b_g)
    dbf_g[1:I_g-1] = b_g[2:I_g] - b_g[1:I_g-1]
    dbf_g[I_g] = dbf_g[I_g-1]

    dbb_g = similar(b_g)
    dbb_g[2:I_g] = b_g[2:I_g] - b_g[1:I_g-1]
    dbb_g[1] = dbb_g[2]

    daf_g = similar(a_g)
    daf_g[1:J_g-1] = a_g[2:J_g] - a_g[1:J_g-1]
    daf_g[J_g] = daf_g[J_g-1]

    dab_g = similar(a_g)
    dab_g[2:J_g] = a_g[2:J_g] - a_g[1:J_g-1]
    dab_g[1] = dab_g[2]

    db_g_tilde       = 0.5*(dbb_g + dbf_g)
    db_g_tilde[1]    = 0.5*dbf_g[1]
    db_g_tilde[end]  = 0.5*dbb_g[end]

    da_g_tilde       = 0.5*(dab_g + daf_g)
    da_g_tilde[1]    = 0.5*daf_g[1]
    da_g_tilde[end]  = 0.5*dab_g[end]

    return kron(da_g_tilde, db_g_tilde)
end

"""
```
@inline function transition(ddeath::R, pam::R, ξ::R, w::S, a_lb::R, aggZ::S,
                            d::Array{U,3}, d_g::Array{U,3}, s::Array{U,3}, s_g::Array{U,3},
                            r_a::U, a::Array{R,1}, a_g::Array{R,1}, b::Array{R,1},
                            b_g::Array{R,1}, y::Vector{U}, cost, util, deposit;
                            permanent::Bool=false) where {R<:AbstractFloat, S<:Number, U<:Number}
```
Constructs matrix which captures transition dynamics for both HJB and KF equations.
"""
@inline function transition(ddeath::R, pam::R, ξ::R, w::S, a_lb::R, aggZ::S, d::Array{U,3},
                            d_g::Array{U,3}, s::Array{U,3}, s_g::Array{U,3}, r_a::U,
                            a::Array{R,1}, a_g::Array{R,1}, b::Array{R,1}, b_g::Array{R,1},
                            y::Vector{U}, cost; permanent::Bool=false) where {R<:AbstractFloat,
                                                                              S<:Number, U<:Number}
    I, J, N     = size(d)
    I_g, J_g, _ = size(d_g)

    perm_const = r_a + (permanent ? ddeath * pam - aggZ : ddeath * pam)
    ξw       = ξ * w

    χ, ζ = similar(d), similar(d)
    X, Z = similar(d), similar(d)

    # Compute drifts for HJB
    for i=1:I, j=1:J, n=1:N

        adrift = a[j] * perm_const + (ξw * y[n])
        bdrift = -d[i,j,n] - cost(d[i,j,n], a[j])

        χ[i,j,n] = (j==1) ? 0.0 : -(min(d[i,j,n], 0) + min(adrift, 0)) / (a[j] - a[j-1])
        ζ[i,j,n] = (j==J) ? 0.0 :  (max(d[i,j,n], 0) + max(adrift, 0)) / (a[j+1] - a[j])

        X[i,j,n] = (i==1) ? 0.0 : -(min(bdrift, 0) + min(s[i,j,n], 0)) / (b[i] - b[i-1])
        Z[i,j,n] = (i==I) ? 0.0 :  (max(bdrift, 0) + max(s[i,j,n], 0)) / (b[i+1] - b[i])
    end

    yyy = -vec(X .+ Z .+ χ .+ ζ)

    χ = reshape(χ, I*J, N)
    χ = circshift(χ, -I)
    ζ = reshape(ζ, I*J, N)
    ζ = circshift(ζ, I)

    X  = reshape(X, I*J, N)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J, N)
    Z  = circshift(Z, 1)

    A = spdiagm(0 => yyy, I => vec(ζ)[I+1:end], -I => vec(χ)[1:end-I],
                          1 => vec(Z)[2:end],   -1 => vec(X)[1:end-1])

    # TODO: is it possible that below (ddeathpam) is incorrect?
    χu, ζu = similar(d_g), similar(d_g)
    Xu, Zu = similar(d_g), similar(d_g)

    # Compute drifts for KFE
    for i=1:I_g, j=1:J_g, n=1:N
        audrift = d_g[i,j,n] + a_g[j] * (r_a + ddeath*pam) + ξw * y[n] -
            (permanent ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j,n] - d_g[i,j,n] - cost(d_g[i,j,n], a_g[j]) -
            (permanent ? aggZ * b_g[i] : 0.0)

        χu[i,j,n] = (j == 1)   ? 0.0 : -min(audrift, 0) / (a_g[j] - a_g[j-1])
        ζu[i,j,n] = (j == J_g) ? 0.0 :  max(audrift, 0) / (a_g[j+1] - a_g[j])

        Xu[i,j,n] = (i == 1)   ? 0.0 : -min(budrift, 0) / (b_g[i] - b_g[i-1])
        Zu[i,j,n] = (i == I_g) ? 0.0 :  max(budrift, 0) / (b_g[i+1] - b_g[i])
    end
    yyyu  = -vec(χu .+ ζu .+ Xu .+ Zu)

    χu = reshape(χu,I_g*J_g,N)
    χu = circshift(χu,-I_g)
    ζu = reshape(ζu,I_g*J_g,N)
    ζu = circshift(ζu,I_g)

    Xu = reshape(Xu, I_g*J_g, N)
    Xu = circshift(Xu, -1)
    Zu = reshape(Zu, I_g*J_g, N)
    Zu = circshift(Zu, 1)

    AT = spdiagm(0 => yyyu, I_g => vec(ζu)[I_g+1:end], -I_g => vec(χu)[1:end-I_g],
                              1 => vec(Zu)[2:end],       -1 => vec(Xu)[1:end-1])
    return A, AT
end

@inline function transition_deriva_lite(permanent::Bool, ddeath::R, pam::R, ξ::R, w::S,
                                        a_lb::R, aggZ::S, d::Array{U,2}, d_g::Array{U,2},
                                        s::Array{U,2}, s_g::Array{U,2}, r_a::U, a::Array{R,1},
                                        a_g::Array{R,1}, b::Array{R,1}, b_g::Array{R,1},
                                        y::U, lambda_ii::R, cost) where {R<:AbstractFloat,
                                                                         S<:Number, U<:Number}
    I, J = size(d)
    I_g, J_g = size(d_g)

    perm_const = permanent ? ddeath * pam - aggZ : ddeath * pam

    χ, ζ = similar(d), similar(d)
    X, Z = similar(d), similar(d)

    for i=1:I, j=1:J
        adrift = a[j] * (r_a + perm_const) + (ξ * w * real(y))

        χ[i,j] = (j==1) ? 0.0 : -(min(d[i,j], 0) + min(adrift, 0)) / (a[j] - a[j-1])
        ζ[i,j] = (j==J) ? 0.0 :  (max(d[i,j], 0) + max(adrift, 0)) / (a[j+1] - a[j])

        X[i,j] = (i==1) ? 0.0 : -(min(-d[i,j] - cost(d[i,j], a[j]), 0) +
                                  min(s[i,j], 0)) / (b[i] - b[i-1])
        Z[i,j] = (i==I) ? 0.0 :  (max(-d[i,j] - cost(d[i,j], a[j]), 0) +
                                  max(s[i,j], 0)) / (b[i+1] - b[i])
    end
    yyy = -vec(X .+ Z .+ χ .+ ζ)

    χ = reshape(χ, I*J)
    χ = circshift(χ, -I)
    ζ = reshape(ζ, I*J)
    ζ = circshift(ζ, I)

    X  = reshape(X, I*J)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J)
    Z  = circshift(Z, 1)

    A = spdiagm(0 => yyy, I => vec(ζ)[I+1:end], -I => vec(χ)[1:end-I],
                          1 => vec(Z)[2:end],   -1 => vec(X)[1:end-1])

    χu, ζu = similar(d_g), similar(d_g)
    Xu, Zu = similar(d_g), similar(d_g)

    for i=1:I_g, j=1:J_g
        # Compute drifts for KFE -- is it possible that below (ddeathpam) is incorrect?
        audrift = d_g[i,j] + a_g[j] * (r_a + ddeath*pam) + ξ * w * y -
            (permanent ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j] - d_g[i,j] - cost(d_g[i,j], max(a_g[j], a_lb)) -
            (permanent ? aggZ * b_g[i] : 0.0)

        χu[i,j] = (j==1)   ? 0.0 : -min(audrift, 0) / (a_g[j] - a_g[j-1])
        ζu[i,j] = (j==J_g) ? 0.0 :  max(audrift, 0) / (a_g[j+1] - a_g[j])

        Xu[i,j] = (i==1)   ? 0.0 : -min(budrift, 0) / (b_g[i] - b_g[i-1])
        Zu[i,j] = (i==I_g) ? 0.0 :  max(budrift, 0) / (b_g[i+1] - b_g[i])
    end
    yyyu  = -vec(χu .+ ζu .+ Xu .+ Zu) .+ lambda_ii

    χu = reshape(χu, I_g*J_g)
    χu = circshift(χu, -I_g)
    ζu = reshape(ζu, I_g*J_g)
    ζu = circshift(ζu, I_g)

    Xu = reshape(Xu, I_g*J_g)
    Xu = circshift(Xu, -1)
    Zu = reshape(Zu, I_g*J_g)
    Zu = circshift(Zu, 1)

    # AT is pre-transposed for efficiency
    AT = spdiagm(0 => yyyu, -I_g => vec(ζu)[I_g+1:end], I_g => vec(χu)[1:end-I_g],
                            -1 => vec(Zu)[2:end],          1 => vec(Xu)[1:end-1])
    return A, AT
end

"""
```
@inline function construct_problem_functions(γ::T, χ0::R, χ1::R, χ2::R,
                                             a_lb::R) where {T<:Number,R<:Number}
```
Initializes utility, adjustment cost, and deposit functions.
"""
@inline function construct_problem_functions(γ::T, χ0::R, χ1::R, χ2::R,
                                             a_lb::R) where {T<:Number,R<:Number}

    util = if γ == 1.0
        @inline util1(c::S) where {S<:Number} = log(c)
    else
        @inline util2(c::S) where {S<:Number} = 1.0 / (1.0 - γ) * (c ^ (1-γ) - 1.0)
    end

    @inline function deposit(Va, Vb, a)
        indx_plus  = ((Va / Vb - 1 - χ0) > 0)
        indx_minus = ((Va / Vb - 1 + χ0) < 0)
        return χ1 * (max(Va / Vb - 1 - χ0, 0)) ^ (1/χ2) * a * indx_plus +
            (-χ1) * (max(-(Va / Vb - 1) - χ0, 0)) ^ (1/χ2) * a * indx_minus
    end

    @inline function cost(d, a)
        d_scaled = abs(d / max(a, a_lb))
        return max(a, a_lb) * (χ0 * (d_scaled) + 1.0 / (1.0 + χ2) *
                               ((d_scaled)^(1+χ2) * χ1^(-χ2)))
    end
    return util, deposit, cost
end

"""
```
@inline function create_y_grid(y_size::Int64, ygrid_new::Int64)
```
Function initializes income grid.
"""
@inline function create_y_grid(y_size::Int64, ygrid_new::Int64)#, dataroot::String)

    y      = Array{Float64}(undef, 0)
    y_dist = Array{Float64}(undef, 0)
    λ      = Array{Float64}(undef, 0)

    if y_size == 2
        y  = [0.1, 1]

        # Transition probabilities
        λ1 = 0.5 # expected duration of unemployment is 2 quarters
        λ2 = (λ1 / (y[2] * .93 - y[1])) * (y[2] - y[2] * .93) # unemployment rate 7
        λ  = [-λ1 λ1; λ2 -λ2]

    elseif y_size == 30
        # Load grid and transition probabilities from external calibration
        #y      = load(dataroot * "income_grid_30.jld2", "y")
        #y_dist = load(dataroot * "income_grid_30.jld2", "y_dist")
        #λ      = load(dataroot * "income_transition_30.jld2", "lambda")
        y = [-4.0043487, -2.3712946, -2.2529091, -1.861575, -1.7594332, -1.751468,
             -1.7514113, -1.7434461, -1.6413042, -1.1315846, -0.61985494, -0.50146947,
             -0.1101354, -0.0079935576, -2.8361989e-5, 2.8362067e-5, 0.0079935649,
             0.11013543, 0.50146947, 0.619855, 1.1315847, 1.6413042, 1.7434461, 1.7514113,
             1.751468, 1.7594332, 1.8615751, 2.2529091, 2.3712946, 4.0043487]

        y_dist = [0.001696604, 0.0036905108, 0.054850618, 0.0034315242, 0.002345745,
                  0.0034002893, 0.0034002799, 0.0023457452, 0.0034315242, 0.0036905106,
                  0.11931293, 0.001696604, 0.11093998, 0.07583712, 0.10993017, 0.10992986,
                  0.075837123, 0.11093998, 0.001696604, 0.11931292, 0.0036905108, 0.0034315242,
                  0.002345745, 0.0034002893, 0.0034002799, 0.0023457452, 0.0034315242,
                  0.054850615, 0.0036905106, 0.001696604]

        λ = build_lambda(y_size)

        # Create grid of income shocks
        y = exp.(y)

    elseif y_size == 33
        # Load grid and transition probabilities from external calibration
        #y      = load(dataroot * "income_grid.jld2", "y")
        #y_dist = load(dataroot * "income_grid.jld2", "y_dist")
        #λ      = load(dataroot * "income_transition.jld2", "lambda")
        y = [-4.7415191, -3.9041665, -3.1769113, -2.8490016, -2.5749613, -2.1240338,
             -2.011649, -1.8925175, -1.6610012, -1.2843938, -1.2100737, -0.95648415,
             -0.68244382, -0.60812364, -0.23151632, -0.1191315, 0.0, 0.1191315, 0.23151632,
             0.60812364, 0.68244382, 0.95648415, 1.2100737, 1.2843938, 1.6610012, 1.8925175,
             2.011649, 2.1240338, 2.5749613, 2.8490016, 3.1769113, 3.9041665, 4.7415191]

        y_dist = [0.00026325521, 0.00073171922, 0.0015341331, 0.008947109, 0.0025562662,
                  0.0032773437, 0.024868536, 0.011062805, 0.0032773437, 0.052139733,
                  0.0025562662, 0.00026325521, 0.086878404, 0.0015341331, 0.11138526,
                  0.00073171922, 0.37598543, 0.00073171922, 0.11138526, 0.0015341331,
                  0.086878404, 0.00026325521, 0.0025562662, 0.052139733, 0.0032773437,
                  0.011062805, 0.024868536, 0.0032773437, 0.0025562662, 0.008947109,
                  0.0015341331, 0.00073171922, 0.00026325521]

        λ = build_lambda(y_size)
        # Create grid of income shocks
        y = exp.(y)
    end

    # Compute stationary income distribution
    y_dist = stat_dist(λ')   # stationary distribution
    y_mean = dot(y, y_dist)

    if     ygrid_new == 0
        y = y ./ y_mean
    elseif ygrid_new == 1
        y = 0.1907858 .* y ./ y_mean
    end

    y_mean = dot(y, y_dist)
    N = size(y, 2)
    return λ, y, y_mean, y_dist, N
end

"""
```
@inline function create_a_grid(agrid_new::Int64, J::Int64, J_g::Int64, amin::Int64, amax::Int64)
```
Function initializes illiquid asset grids.
"""
@inline function create_a_grid(agrid_new::Int64, J::Int64, J_g::Int64, amin::Int64, amax::Int64)
    a   = Array{Float64}(undef, J)
    a_g = Array{Float64}(undef, J_g)

    if agrid_new == 0

        a 		    = range(0, stop = 1, length = J) # set the grid
        coeff_power = 0.9
        power       = 8
        a 			= amax*((1 - coeff_power) * a + coeff_power * (a .^ power))

        a_g         = range(0, stop = 1, length = J_g)
        amax_g      = 150
        coeff_power = 0.9
        power       = 8
        a_g         = amax_g*((1-coeff_power) * a_g + coeff_power * (a_g .^ power))

    elseif agrid_new == 1

        agridparam = 0.15
        a_raw = range(0, stop = 1, length = J)
        a_raw = a_raw.^(1/agridparam)

        a = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a[i] = (i-1) * a[10]/(10 - 1)
        end

        # Define new grid points for g here. For now, it is created the same.
        agridparam = 0.15

        a_raw = range(0, stop = 1, length = J_g)
        a_raw = a_raw.^(1/agridparam)
        a_g = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a_g[i] = (i-1) * a_g[10]/(10 - 1)
        end

    elseif agrid_new == 2

        agridparam = 0.15

        a_raw = range(0, stop = 1, length = J)
        a_raw = a_raw.^(1/agridparam)
        a = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a[i] = (i-1) * a[10]/(10 - 1)
        end

        # Define new grid points for g here. For now, it is created the same.
        agridparam = 0.21

        a_raw = range(0, stop = 1, length = J_g)
        a_raw = a_raw.^(1/agridparam)
        a_g = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a_g[i] = (i-1) * a_g[10]/(10. - 1.)
        end

    end

    a_g_0pos = findall(!iszero, a_g .== 0)

    return a, a_g, a_g_0pos
end

"""
```
@inline function create_b_grid(bgrid_new::Int64, I::Int64)
```
Function to initialize liquid asset grid.
"""
@inline function create_b_grid(bgrid_new::Int64, I::Int64, I_g::Int64)
    I = 50 # TEMP
    I_neg = 10
    I_pos = 40

    bmin = Array{Float64}(undef, 0)
    bmax = Array{Float64}(undef, 0)
    b    = Array{Float64}(undef, 0)
    b_g  = Array{Float64}(undef, 0)

    if bgrid_new == 0

        # Lower bound
        b_    = -2		# lower bound on liquid assets

        # Preparations
        I_neg = Int(I / 10) # number of grid points on the negative part of the real line
        bmin  = b_ 		# lower bound
        bmax  = 80

        # Grid for positive part
        b_pos 			= range(0, stop = 1, length = Int(I-I_neg))
        coeff_power 	= 0.9
        power = 8
        b_pos 			= bmax*((1 - coeff_power) * b_pos + coeff_power * (b_pos.^power))

        # Grid for negative part
        b_neg 			= range(0, stop = 1, length = Int(I_neg+1))
        coeff_power 	= 0
        power           = 1
        b_neg 			= -bmin*((1 - coeff_power) * b_neg + coeff_power * (b_neg.^power)) .+ bmin
        b_neg           = b_neg[1:end-1]

        # Put together grid on liquid assets
        b = [b_neg;b_pos]

        ## Define new grid points for g here. For now, it is created the same.
        I_g_neg = 5
        bmax_g  = 70

        # Grid for positive part
        b_pos 			= range(0, stop = 1, length = Int(I_g-I_g_neg))
        coeff_power 	= 0.9
        power = 8
        b_pos 			= bmax_g*((1 - coeff_power) * b_pos + coeff_power * (b_pos.^power))

        # Grid for negative part
        b_neg 			= range(0, stop = 1, length = Int(I_g_neg+1))
        coeff_power 	= 0
        power = 1
        b_neg = -bmin*((1 - coeff_power) * b_neg + coeff_power * (b_neg.^power)) .+ bmin
        b_neg = b_neg[1:end-1]

        # Put together grid on liquid assets
        b_g = [b_neg;b_pos]

    elseif bgrid_new == 1

        # set parameters
        I_neg = I - I_pos
        bmin  = -1
        b0    = 0
        bmax  = 40
        bgridparam = 0.35
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw .^ (1 / bgridparam)
        bp     = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int64(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn; Vector{Float64}(undef, Int64(I_neg/2-1))]
        for i = Int64(I_neg/2+2) : I_neg
            bn[i] = b0 - (bn[Int64(I_neg)+2-i] - bn[1])
        end

        # Put everything together
        b = [bn; bp]

        ## Define new grid points for g here. For now, it is created the same.
        # set parameters
        I_neg = I_g - I_pos
        bmin  = -1
        b0    = 0
        bmax  = 40
        bgridparam = 0.35
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1/bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int64(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int64(I_neg/2-1))]
        for i = Int64(I_neg/2+2):I_neg
            bn[i] = b0 - (bn[Int64(I_neg+2-i)] - bn[1])
        end

        # Put everything together
        b_g = [bn; bp]

    elseif bgrid_new == 2

        # set parameters
        bmin = -1
        b0   = 0
        bmax = 40
        bgridparam = 0.25
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1/bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int(I_neg/2-1))]

        for i = I_neg/2+2:I_neg
            bn[i] = b0 - (bn[I_neg+2-i] - bn[1])
        end

        # put everything together
        b = [bn;bp]

        ## Define new grid points for g here. For now, it is created the same.
        # set parameters
        I_pos = 40
        I_neg = I_g - I_pos
        bmin = -1
        b0 = 0
        bmax = 40
        bgridparam = 0.6
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1 / bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int(I_neg/2-1))]
        for i = Int(I_neg/2+2):I_neg
            bn[i] = b0 - (bn[Int(I_neg+2-i)] - bn[1])
        end

        # put everything together
        b_g = [bn;bp]
    end
    b_g_0pos = findall(b_g .== 0)

    return I_neg, bmin, bmax, b, b_g, b_g_0pos
end
