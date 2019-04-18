"""
Solves Hamilton-Jacobi-Bellman equation.

"""
@inline function solve_hjb(V::Union{Array{R},Array{S}}, I_g::T, J_g::T, a_lb::R, ggamma::R,
                           permanent::Bool, ddeath::R, pam::R, aggZ::Union{R,S}, xxi::R, tau_I::R,
                           w::Union{R,S}, trans::R, r_b_vec::Union{Vector{R},Vector{S}},
                           y::Vector{U}, a::Vector{R}, b::Vector{R}, cost, util, deposit;
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, T<:Int, U<:Number}
    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------
    I, J, N = size(V)

    perm_c = permanent ? ddeath * pam - aggZ : ddeath * pam
    c0_c   = ((1-xxi) - tau_I) * w

    c, s, d = similar(V), similar(V), similar(V)
    for i=1:I, j=1:J, n=1:N

        c0 = c0_c * real(y[n]) + b[i] * (r_b_vec[i] + perm_c) + trans

        # ---- Liquid Assets, Forward + Backward Difference ----
        VaF = (j==J) ? 0.0 : max((V[i,j+1,n] - V[i,j,n]) / (a[j+1]-a[j]), Vamin)
        VaB = (j==1) ? 0.0 : max((V[i,j,n] - V[i,j-1,n]) / (a[j]-a[j-1]), Vamin)

        # ---- Illiquid Assets, Forward + Backward Difference ----
        VbF = (i==I) ? 0.0 : max((V[i+1,j,n] - V[i,j,n]) / (b[i+1]-b[i]), Vbmin)
        VbB = (i==1) ? 0.0 : max((V[i,j,n] - V[i-1,j,n]) / (b[i]-b[i-1]), Vbmin)

        # Decisions conditional on a particular direction of derivative
        cF = (i==I) ? 0.0 : VbF ^ (-1 / ggamma)
        cB = (i==1) ? c0  : VbB ^ (-1 / ggamma)

        sF = (i==I) ? 0.0 : c0 - cF
        sB = (i==1) ? 0.0 : c0 - cB

        if (i==1 && j > 1) VbB = util(cB) end

        #----------------------------------------------------------------
        # Consumption  & Savings Decision
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

@inline function solve_hjb_lite(V::Union{Array{R},Array{S}}, I_g::T, J_g::T, a_lb::R, ggamma::R,
                           permanent::Bool, ddeath::R, pam::R, aggZ::Union{R,S}, xxi::R, tau_I::R,
                           w::Union{R,S}, trans::R, r_b_vec::Union{Vector{R},Vector{S}},
                           y::U, a::Vector{R}, b::Vector{R}, cost, util, deposit;
                           Vamin = 0.0, Vbmin = 1e-8) where {R<:Float64, S, T<:Int, U<:Number}
    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------
    I, J = size(V)

    perm_c = permanent ? ddeath * pam - aggZ : ddeath * pam
    c0_c   = ((1-xxi) - tau_I) * w

    c, s, d = similar(V), similar(V), similar(V)
    for i=1:I, j=1:J

        c0 = c0_c * real(y) + b[i] * (r_b_vec[i] + perm_c) + trans

        # ---- Liquid Assets, Forward + Backward Difference ----
        VaF = (j==J) ? 0.0 : max((V[i,j+1] - V[i,j]) / (a[j+1]-a[j]), Vamin)
        VaB = (j==1) ? 0.0 : max((V[i,j] - V[i,j-1]) / (a[j]-a[j-1]), Vamin)

        # ---- Illiquid Assets, Forward + Backward Difference ----
        VbF = (i==I) ? 0.0 : max((V[i+1,j] - V[i,j]) / (b[i+1]-b[i]), Vbmin)
        VbB = (i==1) ? 0.0 : max((V[i,j] - V[i-1,j]) / (b[i]-b[i-1]), Vbmin)

        # Decisions conditional on a particular direction of derivative
        cF = (i==I) ? 0.0 : VbF ^ (-1 / ggamma)
        cB = (i==1) ? c0  : VbB ^ (-1 / ggamma)

        sF = (i==I) ? 0.0 : c0 - cF
        sB = (i==1) ? 0.0 : c0 - cB

        if (i==1 && j > 1) VbB = util(cB) end

        #----------------------------------------------------------------
        # Consumption  & Savings Decision
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


@inline function add_diag(A::SparseMatrixCSC, c::AbstractFloat)
    for i=1:size(A,1)
        A[i,i] += c
    end
    return A
end

@inline function update_value_fn(A::SparseMatrixCSC, Delta::R, lambda::Matrix{R},
                                 I::Int64, J::Int64, N::Int64,
                                 u::Union{Array{T,3},Array{R,3}},
                                 Vn::Union{Array{T,3},Array{R,3}},
                                 rrho::R, ddeath::R) where {R<:Float64,T<:Complex}
    Vn_new = Array{Float64,3}(undef,I,J,N)
    for kk = 1:N
        Ak    = A[1+(kk-1)*(I*J):kk*(I*J), 1+(kk-1)*(I*J):kk*(I*J)]
        Bk    = add_diag(-Delta*Ak, 1 + Delta*(rrho + ddeath) - Delta*lambda[kk,kk])
        uk_stacked  = reshape(u[:,:,kk], I * J, 1)
        Vk_stacked  = reshape(Vn[:,:,kk], I * J, 1)
        indx_k      = 1:N .!= kk
        Vkp_stacked = sum(broadcast(*,  lambda[kk, indx_k]',
                                    reshape(Vn[:,:,indx_k], I*J, N-1)), dims=2)
        qk          = Delta .* uk_stacked + Vk_stacked + Delta .* Vkp_stacked
        Vn1k_stacked = Bk \ qk
        Vn_new[:,:,kk]  = reshape(Vn1k_stacked, I, J, 1)
    end
    return Vn_new
end

@inline function set_grids(a, b, a_g, b_g, y, r_b, r_b_borr)
    I   = length(b)
    J   = length(a)
    I_g = length(b_g)
    J_g = length(a_g)
    N   = length(y)

    b_grid    = permutedims(repeat(b, 1, J, N), [1 2 3])
    a_grid    = permutedims(repeat(a, 1, I, N), [2 1 3])
    y_grid    = permutedims(repeat(y, 1, I, J), [2 3 1])
    r_b_grid  = r_b .* (b_grid .>= 0) + r_b_borr .* (b_grid .< 0)

    dbf_grid            = Array{Float64}(undef, I,J,N)
    dbb_grid            = Array{Float64}(undef, I,J,N)
    dbf_grid[1:I-1,:,:] = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbf_grid[I,:,:]     = dbf_grid[I-1,:,:]
    dbb_grid[2:I,:,:]   = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbb_grid[1,:,:]     = dbb_grid[2,:,:]

    daf_grid            = Array{Float64}(undef, I, J, N)
    dab_grid            = Array{Float64}(undef, I, J, N)
    daf_grid[:,1:J-1,:] = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    daf_grid[:,J,:]     = daf_grid[:,J-1,:]
    dab_grid[:,2:J,:]   = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    dab_grid[:,1,:]     = dab_grid[:,2,:]

    db_tilde      = 0.5*(dbb_grid[:,1,1] + dbf_grid[:,1,1])
    db_tilde[1]   = 0.5*dbf_grid[1,1,1]
    db_tilde[end] = 0.5*dbb_grid[end,1,1]
    da_tilde      = 0.5*(dab_grid[1,:,1] + daf_grid[1,:,1])
    da_tilde[1]   = 0.5 * daf_grid[1,1,1]
    da_tilde[end] = 0.5*dab_grid[1,end,1]

    dab_tilde      = kron(da_tilde, db_tilde)
    dab_tilde_grid = reshape(repeat(dab_tilde, N, 1), I, J, N)
    dab_tilde_mat  = spdiagm(0 => vec(repeat(dab_tilde, N, 1)))

    b_g_grid     = permutedims(repeat(b_g, 1, J_g, N),  [1 2 3])
    a_g_grid     = permutedims(repeat(a_g, 1, I_g, N),  [2 1 3])
    y_g_grid     = permutedims(repeat(y, 1, I_g, J_g), [2 3 1])
    r_b_g_grid   = r_b .* (b_g_grid .>= 0) + r_b_borr .* (b_g_grid .< 0)

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

    dab_g_tilde_grid = reshape(repeat(dab_g_tilde,N,1),I_g,J_g,N)

    dab_g_tilde_mat  = spdiagm(0 => vec(repeat(dab_g_tilde,N,1)))

    return a_grid, a_g_grid, b_grid, b_g_grid, y_grid, y_g_grid, r_b_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid
end


@inline function set_grids_lite(a, b, a_g, b_g, N)
    I   = length(b)
    J   = length(a)
    I_g = length(b_g)
    J_g = length(a_g)

    b_grid    = permutedims(repeat(b, 1, J, N), [1 2 3])
    a_grid    = permutedims(repeat(a, 1, I, N), [2 1 3])

    dbf_grid            = Array{Float64}(undef, I,J,N)
    dbb_grid            = Array{Float64}(undef, I,J,N)
    dbf_grid[1:I-1,:,:] = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbf_grid[I,:,:]     = dbf_grid[I-1,:,:]
    dbb_grid[2:I,:,:]   = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbb_grid[1,:,:]     = dbb_grid[2,:,:]

    daf_grid            = Array{Float64}(undef, I, J, N)
    dab_grid            = Array{Float64}(undef, I, J, N)
    daf_grid[:,1:J-1,:] = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    daf_grid[:,J,:]     = daf_grid[:,J-1,:]
    dab_grid[:,2:J,:]   = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    dab_grid[:,1,:]     = dab_grid[:,2,:]

    db_tilde      = 0.5*(dbb_grid[:,1,1] + dbf_grid[:,1,1])
    db_tilde[1]   = 0.5*dbf_grid[1,1,1]
    db_tilde[end] = 0.5*dbb_grid[end,1,1]
    da_tilde      = 0.5*(dab_grid[1,:,1] + daf_grid[1,:,1])
    da_tilde[1]   = 0.5 * daf_grid[1,1,1]
    da_tilde[end] = 0.5*dab_grid[1,end,1]

    dab_tilde      = kron(da_tilde, db_tilde)
    dab_tilde_grid = reshape(repeat(dab_tilde, N, 1), I, J, N)
    dab_tilde_mat  = spdiagm(0 => vec(repeat(dab_tilde, N, 1)))

    b_g_grid     = permutedims(repeat(b_g, 1, J_g, N),  [1 2 3])
    a_g_grid     = permutedims(repeat(a_g, 1, I_g, N),  [2 1 3])

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

    dab_g_tilde_grid = reshape(repeat(dab_g_tilde,N,1),I_g,J_g,N)
    dab_g_tilde_mat  = spdiagm(0 => vec(repeat(dab_g_tilde,N,1)))

    return dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde
end


@inline function set_grids_lite2(a, b, a_g, b_g)
    I   = length(b)
    J   = length(a)
    I_g = length(b_g)
    J_g = length(a_g)

    b_grid    = repeat(b, 1, J)
    a_grid    = permutedims(repeat(a, 1, I), [2 1])

    dbf_grid            = Array{Float64}(undef, I, J)
    dbb_grid            = Array{Float64}(undef, I, J)
    dbf_grid[1:I-1,:] = b_grid[2:I,:] - b_grid[1:I-1,:]
    dbf_grid[I,:]     = dbf_grid[I-1,:]
    dbb_grid[2:I,:]   = b_grid[2:I,:] - b_grid[1:I-1,:]
    dbb_grid[1,:]     = dbb_grid[2,:]

    daf_grid            = Array{Float64}(undef, I, J)
    dab_grid            = Array{Float64}(undef, I, J)
    daf_grid[:,1:J-1] = a_grid[:,2:J] - a_grid[:,1:J-1]
    daf_grid[:,J]     = daf_grid[:,J-1]
    dab_grid[:,2:J]   = a_grid[:,2:J] - a_grid[:,1:J-1]
    dab_grid[:,1]     = dab_grid[:,2]

    db_tilde      = 0.5*(dbb_grid[:,1] + dbf_grid[:,1])
    db_tilde[1]   = 0.5*dbf_grid[1,1]
    db_tilde[end] = 0.5*dbb_grid[end,1]
    da_tilde      = 0.5*(dab_grid[1,:] + daf_grid[1,:])
    da_tilde[1]   = 0.5 * daf_grid[1,1]
    da_tilde[end] = 0.5*dab_grid[1,end]

    dab_tilde      = kron(da_tilde, db_tilde)
    dab_tilde_grid = reshape(dab_tilde, I, J)
    dab_tilde_mat  = spdiagm(0 => vec(dab_tilde))

    b_g_grid     = repeat(b_g, 1, J_g)
    a_g_grid     = permutedims(repeat(a_g, 1, I_g), [2 1])

    dbf_g_grid = Array{Float64}(undef, I_g, J_g)
    dbf_g_grid[1:I_g-1,:] = b_g_grid[2:I_g,:] - b_g_grid[1:I_g-1,:]
    dbf_g_grid[I_g,:] = dbf_g_grid[I_g-1,:]
    dbb_g_grid = Array{Float64}(undef, I_g, J_g)
    dbb_g_grid[2:I_g,:] = b_g_grid[2:I_g,:] - b_g_grid[1:I_g-1,:]
    dbb_g_grid[1,:] = dbb_g_grid[2,:]

    daf_g_grid = Array{Float64}(undef, I_g, J_g)
    daf_g_grid[:,1:J_g-1] = a_g_grid[:,2:J_g] - a_g_grid[:,1:J_g-1]
    daf_g_grid[:,J_g]     = daf_g_grid[:,J_g-1]
    dab_g_grid = Array{Float64}(undef, I_g, J_g)
    dab_g_grid[:,2:J_g] = a_g_grid[:,2:J_g] - a_g_grid[:,1:J_g-1]
    dab_g_grid[:,1]     = dab_g_grid[:,2]

    db_g_tilde       = 0.5*(dbb_g_grid[:,1] + dbf_g_grid[:,1])
    db_g_tilde[1]    = 0.5*dbf_g_grid[1,1]
    db_g_tilde[end]  = 0.5*dbb_g_grid[end,1]
    da_g_tilde       = 0.5*(dab_g_grid[1,:] + daf_g_grid[1,:])
    da_g_tilde[1]    = 0.5*daf_g_grid[1,1]
    da_g_tilde[end]  = 0.5*dab_g_grid[1,end]
    dab_g_tilde      = kron(da_g_tilde, db_g_tilde)

    dab_g_tilde_grid = reshape(dab_g_tilde, I_g, J_g)
    dab_g_tilde_mat  = spdiagm(0 => vec(dab_g_tilde))

    return dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde
end


"""
```
@inline function set_vectors(a, b, a_g, b_g, N, r_b, r_b_borr)
```
Instantiates necessary difference vectors.
"""
@inline function set_vectors(a::Vector{T}, b::Vector{T}, a_g::Vector{T}, b_g::Vector{T},
                             N::Int) where {T<:Float64}

    I, J     = length(b),   length(a)
    I_g, J_g = length(b_g), length(a_g)

    dbf = similar(b)
    dbf[1:I-1] = b[2:I] - b[1:I-1]
    dbf[I]     = dbf[I-1]

    dbb = similar(b)
    dbb[2:I] = b[2:I] - b[1:I-1]
    dbb[1]   = dbb[2]

    daf = similar(a)
    daf[1:J-1] = a[2:J] - a[1:J-1]
    daf[J]     = daf[J-1]

    dab = similar(a)
    dab[2:J] = a[2:J] - a[1:J-1]
    dab[1]   = dab[2]

    db_tilde = 0.5*(dbb + dbf)
    db_tilde[1] = 0.5*(dbf[1])
    db_tilde[end] = 0.5*(dbb[end])

    da_tilde = 0.5*(dab + daf)
    da_tilde[1] = 0.5*(daf[1])
    da_tilde[end] = 0.5*(dab[end])

    dab_tilde      = kron(da_tilde, db_tilde)
    dab_tilde_grid = reshape(repeat(dab_tilde, N, 1), I, J, N)
    dab_tilde_mat  = spdiagm(0 => vec(repeat(dab_tilde, N, 1)))

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

    dab_g_tilde      = kron(da_g_tilde, db_g_tilde)
    dab_g_tilde_grid = reshape(repeat(dab_g_tilde,N,1),I_g,J_g,N)
    dab_g_tilde_mat  = spdiagm(0 => vec(repeat(dab_g_tilde,N,1)))

    return dab, dab_g, dab_tilde, dab_g_tilde, dab_tilde_grid, dab_tilde_mat, dab_g_tilde_grid, dab_g_tilde_mat
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

    dab_g_tilde      = kron(da_g_tilde, db_g_tilde)

    return dab_g_tilde
end


@inline function transition(ddeath::S, pam::S, xxi::S, w::R, chi0::S, chi1::S, chi2::S, a_lb::S,
                            r_a::R, y::Vector{T}, d::Union{Array{S,3},Array{T,3}},
                            d_g::Union{Array{S,3},Array{T,3}},
                            s::Union{Array{S,3},Array{T,3}}, s_g::Union{Array{S,3},Array{T,3}},
                            a::Vector{S}, a_g::Vector{S}, b::Vector{S},
                            b_g::Vector{S}) where {S<:AbstractFloat, R<:Number, T<:Number}
    I, J, N  = size(d)
    I_g, J_g = size(d_g)

    chi  = Array{Float64,3}(undef, I,J,N)
    yy   = Array{Float64,3}(undef, I,J,N)
    zeta = Array{Float64,3}(undef, I,J,N)

    X = Array{Float64,3}(undef, I,J,N)
    Y = Array{Float64,3}(undef, I,J,N)
    Z = Array{Float64,3}(undef, I,J,N)

    chiu  = Array{Float64,3}(undef, I_g,J_g,N)
    yyu   = Array{Float64,3}(undef, I_g,J_g,N)
    zetau = Array{Float64,3}(undef, I_g,J_g,N)

    Xu = Array{Float64,3}(undef, I_g,J_g,N)
    Yu = Array{Float64,3}(undef, I_g,J_g,N)
    Zu = Array{Float64,3}(undef, I_g,J_g,N)

    α    = (r_a + ddeath*pam)
    xxiw = xxi * w

    for i=1:I, j=1:J, n=1:N
        # compute all drifts
        adrift = a[j] * α + xxiw * y[n]
        bdrift = -d[i,j,n] - adj_cost_fn(d[i,j,n], a[j], chi0, chi1, chi2, a_lb)

        chi[i,j,n]  = (j==1) ? 0.0 : -(min(d[i,j,n],0) + min(adrift, 0)) / (a[j] - a[j-1])
        zeta[i,j,n] = (j==J) ? 0.0 :  (max(d[i,j,n],0) + max(adrift, 0)) / (a[j+1] - a[j])

        X[i,j,n] = (i==1) ? 0.0 : -(min(bdrift, 0) + min(s[i,j,n], 0)) / (b[i] - b[i-1])
        Z[i,j,n] = (i==I) ? 0.0 :  (max(bdrift, 0) + max(s[i,j,n], 0)) / (b[i+1] - b[i])
    end
    yy = -chi .- zeta
    Y  = -X .- Z

    centdiag = reshape(yy,I*J,N)
    lowdiag  = reshape(chi,I*J,N)
    lowdiag  = circshift(lowdiag,-I)
    updiag   = reshape(zeta,I*J,N)
    updiag   = circshift(updiag,I)

    centdiag = reshape(centdiag,I*J*N,1)
    updiag   = reshape(updiag,I*J*N,1)
    lowdiag  = reshape(lowdiag,I*J*N,1)

    aa = spdiagm(-I => vec(lowdiag)[1:end-I], 0 => vec(centdiag), I => vec(updiag)[I+1:end])

    centdiag = reshape(Y,I*J,N)
    lowdiag  = reshape(X,I*J,N)
    lowdiag  = circshift(lowdiag,-1)
    updiag   = reshape(Z,I*J,N)
    updiag   = circshift(updiag,1)

    centdiag = reshape(centdiag,I*J*N,1)
    updiag   = reshape(updiag,I*J*N,1)
    lowdiag  = reshape(lowdiag,I*J*N,1)

    bb = spdiagm(0 => vec(centdiag), 1 => vec(updiag)[2:end], -1 => vec(lowdiag)[1:end-1])

    for i=1:I_g, j=1:J_g, n=1:N
        audrift = d_g[i,j,n] + a_g[j] * α + xxiw * y[n]
        budrift = s_g[i,j,n] - d_g[i,j,n] -
            adj_cost_fn(d_g[i,j,n], a_g[j], chi0, chi1, chi2, a_lb)

        chiu[i,j,n]  = (j==1)   ? 0.0 : -min(audrift,0) / (a_g[j] - a_g[j-1])
        zetau[i,j,n] = (j==J_g) ? 0.0 :  max(audrift,0) / (a_g[j+1] - a_g[j])

        Xu[i,j,n] = (i==1)      ? 0.0 : -min(budrift,0) / (b_g[i] - b_g[i-1])
        Zu[i,j,n] = (i==I_g)    ? 0.0 :  max(budrift,0) / (b_g[i+1] - b_g[i])
    end

    yyu = -(chiu .+ zetau)
    Yu  = -(Xu .+ Zu)

    centdiagu = reshape(yyu, I_g * J_g, N)
    lowdiagu  = reshape(chiu,I_g * J_g, N)
    lowdiagu  = circshift(lowdiagu, -I_g)
    updiagu   = reshape(zetau, I_g * J_g, N)
    updiagu   = circshift(updiagu, I_g)

    centdiagu = reshape(centdiagu,I_g * J_g * N, 1)
    updiagu   = reshape(updiagu,  I_g * J_g * N, 1)
    lowdiagu  = reshape(lowdiagu, I_g * J_g * N, 1)

    aau = spdiagm(0 => vec(centdiagu), -I_g => vec(lowdiagu)[1:end-I_g],
                  I_g => vec(updiagu)[I_g+1:end])

    centdiagu = reshape(Yu,I_g*J_g,N)
    lowdiagu = reshape(Xu,I_g*J_g,N)
    lowdiagu = circshift(lowdiagu,-1)
    updiagu = reshape(Zu,I_g*J_g,N)
    updiagu = circshift(updiagu,1)

    centdiagu = reshape(centdiagu,I_g*J_g*N,1)
    updiagu   = reshape(updiagu,I_g*J_g*N,1)
    lowdiagu  = reshape(lowdiagu,I_g*J_g*N,1)

    bbu = spdiagm(0 => vec(centdiagu), 1 => vec(updiagu)[2:end], -1 => vec(lowdiagu)[1:end-1])

    return aa, bb, aau, bbu
end

@inline function transition_deriva2(perm::T, ddeath::R, pam::R, xxi::R, w::S,
                                    a_lb::R, d, d_g, s, s_g, r_a, aggZ::S,
                                    a, a_g, b, b_g, y, lambda,
                                    cost, util, deposit) where {R<:AbstractFloat, S, T}
    I, J, N  = size(d)
    I_g, J_g = size(d_g)
    permanent = (perm==1) ? true : false

    # Compute drifts for HJB
    perm_const = permanent ? ddeath * pam - aggZ : ddeath * pam

    X, Z = similar(d), similar(d)
    #aa   = spzeros(eltype(d), I*J*N, I*J*N)
    #bb   = spzeros(eltype(d), I*J*N, I*J*N)
    #A   = kron(lambda, my_speye(eltype(d), I*J))
    A = #=kron(lambda, my_speye(I*J))=#spzeros(eltype(d), I*J*N, I*J*N)
    f_ind(x, my_J, n) = (x + (n-1)) % my_J + 1
    #chi_dict  = Dict(f_ind.(1:J, J, 1)   .=> 1:J)
    #zeta_dict = Dict(f_ind.(1:J, J, J-1) .=> 1:J)

    #for n=1:N
    #    A[1+(n-1)*I*J:n*I*J, 1+(n-1)*I*J:n*I*J] = lambda
    #end

    for i=1:I, j=1:J, n=1:N

        α = a[j] * (r_a + perm_const) + (xxi * w * real(y[n]))

        chi  = (j==1) ? 0.0 : -(min(norm(d[i,j,n]), 0) + min(α, 0)) / (a[j]   - a[j-1])
        zeta = (j==J) ? 0.0 :  (max(norm(d[i,j,n]), 0) + max(α, 0)) / (a[j+1] - a[j])

        ind      = (I*J)*(n-1) + I*(j-1) + i
        chi_ind  = I*J*(n-1) + I*(f_ind(j,J,J-1) - 1) + i#(chi_dict[j]-1) + i
        zeta_ind = I*J*(n-1) + I*(f_ind(j,J,1) - 1) + i#(zeta_dict[j]-1) + i

        A[ind, ind] += -(chi + zeta)

        if (chi_ind  <= I*J*N - I) A[chi_ind  + I, chi_ind]  += chi  end
        if (zeta_ind >= I + 1)     A[zeta_ind - I, zeta_ind] += zeta end

        #X[i,j,n] = (i==1) ? 0.0 : -(min(-d[i,j,n] -
        #                          adj_cost_fn(d[i,j,n], a[j], chi0, chi1, chi2, a_lb), 0) +
        #                      min(s[i,j,n], 0)) / (b[i] - b[i-1])
        #Z[i,j,n] = (i==I) ? 0.0 : (max(-d[i,j,n] -
        #                         adj_cost_fn(d[i,j,n], a[j], chi0, chi1, chi2, a_lb), 0) +
        #                     max(s[i,j,n], 0)) / (b[i+1] - b[i])

        X = (i==1) ? 0.0 : -(min(-d[i,j,n] - cost(d[i,j,n], max(a[j], a_lb)), 0) +
                              min(s[i,j,n], 0)) / (b[i] - b[i-1])
        Z = (i==I) ? 0.0 : (max(-d[i,j,n] - cost(d[i,j,n], max(a[j], a_lb)), 0) +
                             max(s[i,j,n], 0)) / (b[i+1] - b[i])

        A[ind, ind] += -(X + Z)

        ind      = (I*J)*(n-1) + I*(j-1) + i
        X_ind = I*J*(n-1) + f_ind(I*(j-1) + i, I*J,I*J-1) #I*(chi_dict[j]-1) + i
        Z_ind = I*J*(n-1) + f_ind(I*(j-1) + i, I*J,1)     #I*(zeta_dict[j]-1) + i
        if (X_ind <= I*J*N - 1) A[X_ind + 1, X_ind] += X end
        if (Z_ind >= 2)         A[Z_ind - 1, Z_ind] += Z end
    end
#=
    Y = -(X .+ Z)

    X  = reshape(X, I*J, N)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J, N)
    Z  = circshift(Z, 1)

    bb = spdiagm(0 => vec(Y), 1 => vec(Z)[2:end], -1 => vec(X)[1:end-1])
    @assert bb == bb2
=#
    #aau = spzeros(eltype(d), I_g*J_g*N, I_g*J_g*N)
    #bbu = spzeros(eltype(d), I_g*J_g*N, I_g*J_g*N)
    AT = spzeros(eltype(d), I_g*J_g*N, I_g*J_g*N)

    for i=1:I_g, j=1:J_g, n=1:N

        # Compute drifts for KFE -- is it possible that below (ddeathpam) is incorrect?
        audrift = d_g[i,j,n] + a_g[j] * (r_a + ddeath*pam) + xxi * w * y[n] -
            (permanent ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j,n] - d_g[i,j,n] - cost(d_g[i,j,n], max(a_g[j], a_lb)) -
            (permanent ? aggZ * b_g[i] : 0.0)

        chiu  = (j==1)   ? 0.0 : -min(audrift, 0) / (a_g[j]   - a_g[j-1])
        zetau = (j==J_g) ? 0.0 :  max(audrift, 0) / (a_g[j+1] - a_g[j])

        ind       = (I_g*J_g)*(n-1) + I_g*(j-1) + i
        chiu_ind  = I_g*J_g*(n-1) + I_g *(f_ind(j,J_g,J_g-1)-1) + i
        zetau_ind = I_g*J_g*(n-1) + I_g*(f_ind(j,J_g,1)-1) + i

        AT[ind, ind] = -(chiu + zetau)
        if (chiu_ind  <= I_g*J_g*N - I_g) AT[chiu_ind  + I_g, chiu_ind]  = chiu  end
        if (zetau_ind >= I_g + 1)         AT[zetau_ind - I_g, zetau_ind] = zetau end

        Xu = (i==1)   ? 0.0 : -min(budrift, 0) / (b_g[i]   - b_g[i-1])
        Zu = (i==I_g) ? 0.0 :  max(budrift, 0) / (b_g[i+1] - b_g[i])

        Xu_ind = I_g*J_g*(n-1) + f_ind(I_g*(j-1) + i, I_g*J_g, I_g*J_g - 1)
        Zu_ind = I_g*J_g*(n-1) + f_ind(I_g*(j-1) + i, I_g*J_g, 1)

        AT[ind, ind] += -(Xu + Zu)
        if (Xu_ind <= I_g*J_g*N - 1) AT[Xu_ind + 1, Xu_ind] += Xu  end
        if (Zu_ind >= 2)             AT[Zu_ind - 1, Zu_ind] += Zu end

    end
    return A, AT
end

@inline function transition_deriva(permanent::T, ddeath::R, pam::R, xxi::R, w::S,
                                   a_lb::R, aggZ::S, d::Array{U,3}, d_g::Array{U,3},
                                   s::Array{U,3}, s_g::Array{U,3},
                                   r_a::U, a::Array{R,1}, a_g::Array{R,1}, b::Array{R,1},
                                   b_g::Array{R,1}, y::Vector{U},
                                   cost, util, deposit) where {R<:AbstractFloat, S<:Number,
                                                               T<:Bool, U<:Number}
    I, J, N     = size(d)
    I_g, J_g, _ = size(d_g)

    perm_const = permanent ? ddeath * pam - aggZ : ddeath * pam

    chi, zeta = similar(d), similar(d)
    X, Z      = similar(d), similar(d)

    for i=1:I, j=1:J, n=1:N
        α = a[j] * (r_a + perm_const) + (xxi * w * real(y[n]))

        chi[i,j,n]  = (j==1) ? 0.0 : -(min(d[i,j,n], 0) + min(α, 0)) / (a[j]   - a[j-1])
        zeta[i,j,n] = (j==J) ? 0.0 :  (max(d[i,j,n], 0) + max(α, 0)) / (a[j+1] - a[j])

        X[i,j,n] = (i==1) ? 0.0 : -(min(-d[i,j,n] - cost(d[i,j,n], max(a[j], a_lb)), 0) +
                                    min(s[i,j,n], 0)) / (b[i] - b[i-1])
        Z[i,j,n] = (i==I) ? 0.0 :  (max(-d[i,j,n] - cost(d[i,j,n], max(a[j], a_lb)), 0) +
                                    max(s[i,j,n], 0)) / (b[i+1] - b[i])
    end
    yyy = -vec(X .+ Z .+ chi .+ zeta)

    chi  = reshape(chi, I*J, N)
    chi  = circshift(chi, -I)
    zeta = reshape(zeta, I*J, N)
    zeta = circshift(zeta, I)

    X  = reshape(X, I*J, N)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J, N)
    Z  = circshift(Z, 1)

    A = spdiagm(0 => yyy, I => vec(zeta)[I+1:end], -I => vec(chi)[1:end-I],
                          1 => vec(Z)[2:end],      -1 => vec(X)[1:end-1])

    chiu, zetau = similar(d_g), similar(d_g)
    Xu, Zu      = similar(d_g), similar(d_g)
    for i=1:I_g, j=1:J_g, n=1:N
        # Compute drifts for KFE -- is it possible that below (ddeathpam) is incorrect?
        audrift = d_g[i,j,n] + a_g[j] * (r_a + ddeath*pam) + xxi * w * y[n] -
            (permanent ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j,n] - d_g[i,j,n] - cost(d_g[i,j,n], max(a_g[j], a_lb)) -
            (permanent ? aggZ * b_g[i] : 0.0)

        chiu[i,j,n]  = (j==1)   ? 0.0 : -min(audrift, 0) / (a_g[j]   - a_g[j-1])
        zetau[i,j,n] = (j==J_g) ? 0.0 :  max(audrift, 0) / (a_g[j+1] - a_g[j])

        Xu[i,j,n] = (i==1)   ? 0.0 : -min(budrift, 0) / (b_g[i]   - b_g[i-1])
        Zu[i,j,n] = (i==I_g) ? 0.0 :  max(budrift, 0) / (b_g[i+1] - b_g[i])
    end
    yyyu  = -vec(chiu .+ zetau .+ Xu .+ Zu)

    chiu  = reshape(chiu,I_g*J_g,N)
    chiu  = circshift(chiu,-I_g)
    zetau = reshape(zetau,I_g*J_g,N)
    zetau = circshift(zetau,I_g)

    Xu = reshape(Xu, I_g*J_g, N)
    Xu = circshift(Xu, -1)
    Zu = reshape(Zu, I_g*J_g, N)
    Zu = circshift(Zu, 1)

    AT = spdiagm(0 => yyyu, I_g => vec(zetau)[I_g+1:end], -I_g => vec(chiu)[1:end-I_g],
                            1 => vec(Zu)[2:end],          -1 => vec(Xu)[1:end-1])
    return A, AT
end

@inline function transition_deriva_lite(permanent::Bool, ddeath::R, pam::R, xxi::R, w::S,
                                        a_lb::R, aggZ::S, d::Array{U,2}, d_g::Array{U,2},
                                        s::Array{U,2}, s_g::Array{U,2},
                                        r_a::U, a::Array{R,1}, a_g::Array{R,1}, b::Array{R,1},
                                        b_g::Array{R,1}, y::U, lambda_ii::R,
                                        cost, util, deposit) where {R<:AbstractFloat, S<:Number,
                                                                    U<:Number}
    I, J = size(d)
    I_g, J_g = size(d_g)

    perm_const = permanent ? ddeath * pam - aggZ : ddeath * pam

    chi, zeta = similar(d), similar(d)
    X, Z      = similar(d), similar(d)
    for i=1:I, j=1:J
        α = a[j] * (r_a + perm_const) + (xxi * w * real(y))

        chi[i,j]  = (j==1) ? 0.0 : -(min(d[i,j], 0) + min(α, 0)) / (a[j] - a[j-1])
        zeta[i,j] = (j==J) ? 0.0 :  (max(d[i,j], 0) + max(α, 0)) / (a[j+1] - a[j])

        X[i,j] = (i==1) ? 0.0 : -(min(-d[i,j] - cost(d[i,j], max(a[j], a_lb)), 0) +
                                    min(s[i,j], 0)) / (b[i] - b[i-1])
        Z[i,j] = (i==I) ? 0.0 :  (max(-d[i,j] - cost(d[i,j], max(a[j], a_lb)), 0) +
                                    max(s[i,j], 0)) / (b[i+1] - b[i])
    end
    yyy = -vec(X .+ Z .+ chi .+ zeta)

    chi  = reshape(chi, I*J)
    chi  = circshift(chi, -I)
    zeta = reshape(zeta, I*J)
    zeta = circshift(zeta, I)

    X  = reshape(X, I*J)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J)
    Z  = circshift(Z, 1)

    A = spdiagm(0 => yyy, I => vec(zeta)[I+1:end], -I => vec(chi)[1:end-I],
                          1 => vec(Z)[2:end],      -1 => vec(X)[1:end-1])

    chiu, zetau = similar(d_g), similar(d_g)
    Xu, Zu      = similar(d_g), similar(d_g)
    for i=1:I_g, j=1:J_g
        # Compute drifts for KFE -- is it possible that below (ddeathpam) is incorrect?
        audrift = d_g[i,j] + a_g[j] * (r_a + ddeath*pam) + xxi * w * y -
            (permanent ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j] - d_g[i,j] - cost(d_g[i,j], max(a_g[j], a_lb)) -
            (permanent ? aggZ * b_g[i] : 0.0)

        chiu[i,j]  = (j==1)   ? 0.0 : -min(audrift, 0) / (a_g[j] - a_g[j-1])
        zetau[i,j] = (j==J_g) ? 0.0 :  max(audrift, 0) / (a_g[j+1] - a_g[j])

        Xu[i,j] = (i==1)   ? 0.0 : -min(budrift, 0) / (b_g[i] - b_g[i-1])
        Zu[i,j] = (i==I_g) ? 0.0 :  max(budrift, 0) / (b_g[i+1] - b_g[i])
    end
    yyyu  = -vec(chiu .+ zetau .+ Xu .+ Zu) .+ lambda_ii

    chiu  = reshape(chiu,I_g*J_g)
    chiu  = circshift(chiu,-I_g)
    zetau = reshape(zetau,I_g*J_g)
    zetau = circshift(zetau,I_g)

    Xu = reshape(Xu, I_g*J_g)
    Xu = circshift(Xu, -1)
    Zu = reshape(Zu, I_g*J_g)
    Zu = circshift(Zu, 1)

    AT = spdiagm(0 => yyyu, -I_g => vec(zetau)[I_g+1:end], I_g => vec(chiu)[1:end-I_g],
                            -1 => vec(Zu)[2:end],          1 => vec(Xu)[1:end-1])
    return A, AT
end


@inline function construct_problem_functions(γ::T, χ0::R, χ1::R, χ2::R,
                                             a_lb::R) where {T<:Number,R<:Number}

    @inline util(c::S) where {S<:Number} = 1.0 / (1.0 - γ) * (c ^ (1-γ) - 1.0)
    if γ == 1.0
        @inline util(c::S) where {S<:Number} = log(c)
    end

    @inline function deposit(Va, Vb, a)
        indx_plus  = ((Va / Vb - 1 - χ0) > 0)
        indx_minus = ((Va / Vb - 1 + χ0) < 0)
        return χ1 * (max(Va / Vb - 1 - χ0, 0)) ^ (1/χ2) * a * indx_plus +
            (-χ1) * (max(-(Va / Vb - 1) - χ0, 0)) ^ (1/χ2) * a * indx_minus
    end

    @inline function cost(d, a)
        d_scaled = abs(d / max(a, a_lb))
        return max(a, a_lb) * (χ0 * (d_scaled) + 1.0 / (1.0 + χ2) * ((d_scaled)^(1+χ2) * χ1^(-χ2)))
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
    a   = Array{Float64}(undef, 0)
    a_g = Array{Float64}(undef, 0)

    if agrid_new == 0

        a 		    = range(0, stop = 1, length = J) # set the grid
        coeff_power = 0.9
        power       = 8
        a 			= amax*((1 - coeff_power) * a + coeff_power * (a.^power))

        a_g         = range(0, stop = 1, length = J_g)
        amax_g      = 150
        coeff_power = 0.9
        power       = 8
        a_g         = amax_g*((1-coeff_power) * a_g + coeff_power * (a_g.^power))

    elseif agrid_new == 1

        agridparam = 0.15
        a_raw = range(0, stop = 1, length = J)
        a_raw = a_raw.^(1/agridparam)

        a = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a[i] = (i-1) * a[10]/(10 - 1)
        end

        ## Define new grid points for g here. For now, it is created the same.
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

        ## Define new grid points for g here. For now, it is created the same.
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

# This file contains additional helper functions for computing the steady state.

"""
```
@inline function construct_asset_grid(I::Int64, grid_param::Int64, grid_min::Float64,
                                      grid_max::Float64)
```
Discretize a state space dimension of I gridpoints, where agridparam is the "bending"
coefficient of the grid; i.e. grid_param = 1 implies a uniform grid, and grid_min/grid_max
are the start and endpoints of the grid.
"""
@inline function construct_asset_grid(I::Int64, grid_param::Int64, grid_min::Float64,
                                      grid_max::Float64)
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
